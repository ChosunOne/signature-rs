//! Compiled evaluation plan for the BCH concatenation expression: shared
//! subtrees interned to node ids, evaluated as a straight-line loop over
//! reused per-node buffers instead of recursive per-fold evaluation.

use std::collections::HashMap;
use std::fmt;
use std::hash::Hash;
use std::ops::{AddAssign, MulAssign, Neg};
use std::sync::{Arc, OnceLock};

use commutator_rs::CommutatorTerm;
use lie_rs::{ClassOrder, ClassOrderedCommutation, GatingCache, KernelJob, LieSeries};
use lyndon_rs::generators::Generator;
use num_traits::{One, Zero};

/// One node of the compiled plan.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum DagNode {
    /// Leaf: `0` is the accumulator (`A`), `1` the displacement (`B`).
    Atom(u8),
    /// Internal bracket node; children have strictly smaller ids.
    Binary { left: u32, right: u32 },
}

/// Where a term's coefficient slice lives after evaluation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum TermSource {
    /// The term's root node result, held in the DAG's buffers.
    Node(u32),
    /// The term is the bare displacement atom `B`; its coefficients are the
    /// fold's `rhs` input slice.
    Displacement,
}

pub(crate) struct DagStructure<U> {
    /// Topologically sorted: ids 0 and 1 are the atoms, every internal node's
    /// children have strictly smaller ids.
    pub(crate) nodes: Vec<DagNode>,
    /// Nodes grouped by dependency height (atoms at height 0, a bracket one
    /// above its tallest child). All nodes of a level depend only on earlier
    /// levels, so a level evaluates as one conflict-free parallel batch.
    pub(crate) levels: Vec<Vec<u32>>,
    /// `(source, BCH weight)` per accumulated term, in the original term
    /// order (float summation order preserved).
    terms: Vec<(TermSource, U)>,
}

/// Compiled evaluation plan for one `BchSeriesGenerator` configuration.
///
/// Constructed once per `LogSignatureBuilder::build`; shared across
/// `LogSignature` clones. Scratch buffers are neither cloned nor shared —
/// they are (re)allocated on the first fold after construction or clone.
pub(crate) struct CommutatorDag<U> {
    pub(crate) structure: Arc<DagStructure<U>>,
    /// Result buffer per internal node (node id `i` ↔ buffer `i - 2`).
    pub(crate) buffers: Vec<Vec<U>>,
    /// Scatter-target superset per internal node, mirrors `buffers`.
    ///
    /// A node's list is exactly the set of indices its kernel call scattered
    /// onto when it was last (re)built. That set depends only on the
    /// children's lists — it is a fixed point of the DAG's support
    /// propagation — so the whole collection stays valid while the atom
    /// supports are unchanged, and every list remains a superset of the
    /// node buffer's non-zero indices (the kernel's presence tests only
    /// need the superset: a listed zero-valued index contributes zero).
    pub(crate) nonzeros: Vec<Vec<usize>>,
    /// Deduplication scratch for the collecting rebuild, one bitset per
    /// internal node sized to the basis.
    dirty: Vec<Vec<u64>>,
    /// The atom supports the current `nonzeros` were built from.
    atom_a: Vec<usize>,
    atom_b: Vec<usize>,
    /// Whether `nonzeros` holds a built (or inherited) fixed point.
    lists_built: bool,
    /// Memoized kernel gating keyed by degree-support masks. Valid because
    /// the node lists are a fixed point: while the atom supports are
    /// unchanged, every node's support stays in the same degree slice, so
    /// the (a-mask, b-mask) keys — and hence the active-unit lists — repeat
    /// fold after fold. Built for this DAG's own series/table only.
    gating_cache: GatingCache,
    /// The basis' class-contiguous ordering, created on the first fold and
    /// shared by every kernel call of every fold through this DAG (and by
    /// `clone_shallow` copies): the O(basis) planning is paid once.
    ///
    /// All fold work runs in this ordering — node buffers are class-ordered
    /// and the support lists class-indexed — and
    /// [`Self::accumulate_terms`] applies the single public-order epilogue.
    class_order: OnceLock<Arc<ClassOrder>>,
}

/// Coefficient slice for node `idx`: an atom reads the fold's inputs, an
/// internal node its own result buffer (which precedes `before`'s end since
/// children are topologically earlier).
#[inline]
fn node_slice<'a, U>(idx: u32, before: &'a [Vec<U>], a: &'a [U], b: &'a [U]) -> &'a [U] {
    if idx < 2 {
        if idx == 0 { a } else { b }
    } else {
        &before[idx as usize - 2]
    }
}

/// Non-zero index list counterpart of [`node_slice`].
#[inline]
fn node_nonzeros<'a>(
    idx: u32,
    before: &'a [Vec<usize>],
    a_nonzero: &'a [usize],
    b_nonzero: &'a [usize],
) -> &'a [usize] {
    if idx < 2 {
        if idx == 0 { a_nonzero } else { b_nonzero }
    } else {
        &before[idx as usize - 2]
    }
}

impl<U> CommutatorDag<U>
where
    U: Clone + Zero,
{
    /// Compiles the plan from a BCH series' bracket trees.
    ///
    /// Only terms with non-zero BCH weights are kept; their relative order
    /// (and therefore the accumulation's float summation order) is preserved.
    pub(crate) fn from_bch_series(bch_series: &LieSeries<u8, U>) -> Self
    where
        U: Eq + Hash + One,
    {
        let mut nodes = vec![DagNode::Atom(0), DagNode::Atom(1)];
        // Intern map keyed by the structural `unit_hash`; the bucket vector
        // resolves (astronomically rare) hash collisions by structural
        // equality, so interning never clones a tree for a map key.
        let mut interned: HashMap<u64, Vec<(CommutatorTerm<U, u8>, u32)>> = HashMap::new();

        fn intern<U: Eq + Hash + Clone + One>(
            term: &CommutatorTerm<U, u8>,
            nodes: &mut Vec<DagNode>,
            interned: &mut HashMap<u64, Vec<(CommutatorTerm<U, u8>, u32)>>,
        ) -> u32 {
            match term {
                CommutatorTerm::Atom { atom, .. } => *atom as u32,
                CommutatorTerm::Expression { .. } => {
                    let hash = term.unit_hash();
                    if let Some(bucket) = interned.get(&hash) {
                        for (existing, id) in bucket {
                            if existing == term {
                                return *id;
                            }
                        }
                    }
                    let left =
                        intern(term.left().expect("expression has a left child"), nodes, interned);
                    let right = intern(
                        term.right().expect("expression has a right child"),
                        nodes,
                        interned,
                    );
                    let id = nodes.len() as u32;
                    nodes.push(DagNode::Binary { left, right });
                    interned.entry(hash).or_default().push((term.clone(), id));
                    id
                }
            }
        }

        let mut terms = Vec::new();
        for (i, term) in bch_series.commutator_basis.iter().enumerate() {
            if i == 0 {
                continue; // the degree-1 `A` term: its contribution is the accumulator
            }
            let weight = &bch_series.coefficients[i];
            if weight.is_zero() {
                continue;
            }
            let root = intern(term, &mut nodes, &mut interned);
            let source = match root {
                0 => continue, // defensive: an atom-`A` root adds nothing new
                1 => TermSource::Displacement,
                id => TermSource::Node(id),
            };
            terms.push((source, weight.clone()));
        }

        // Dependency heights and levels: atoms at height 0, a bracket one
        // above its tallest child. The nodes are topologically sorted, so a
        // single forward pass suffices.
        let mut heights = vec![0u32; nodes.len()];
        let mut max_height = 0u32;
        for (idx, node) in nodes.iter().enumerate() {
            if let DagNode::Binary { left, right } = node {
                let h = 1 + heights[*left as usize].max(heights[*right as usize]);
                heights[idx] = h;
                max_height = max_height.max(h);
            }
        }
        let mut levels: Vec<Vec<u32>> = vec![Vec::new(); max_height as usize + 1];
        for (idx, h) in heights.iter().enumerate() {
            levels[*h as usize].push(idx as u32);
        }

        Self {
            structure: Arc::new(DagStructure { nodes, terms, levels }),
            buffers: Vec::new(),
            nonzeros: Vec::new(),
            dirty: Vec::new(),
            atom_a: Vec::new(),
            atom_b: Vec::new(),
            lists_built: false,
            gating_cache: GatingCache::default(),
            class_order: OnceLock::new(),
        }
    }
}

impl<U> CommutatorDag<U>
where
    U: Clone + Default + One + Zero + Eq + MulAssign + Neg<Output = U> + Hash + AddAssign,
{
    /// Evaluates the whole plan for one concatenation.
    ///
    /// `a` is the accumulated log-signature's coefficient slice *before* this
    /// fold's BCH additions, `b` the segment displacement; after the call the
    /// per-term results are available via [`Self::terms`]/[`Self::buffer`]
    /// for accumulation into the log signature.
    pub(crate) fn evaluate<T>(
        &mut self,
        series: &LieSeries<T, U>,
        a: &[U],
        a_nonzero: &[usize],
        b: &[U],
        b_nonzero: &[usize],
    ) where
        T: Clone + Ord + Generator + Hash + Eq,
        U: Send + Sync,
    {
        let internal = self.structure.nodes.len() - 2;
        self.ensure_buffers(series.coefficients.len(), internal);

        // The class-contiguous ordering: created on the first fold, shared
        // by every kernel call of every fold through this DAG. The fold's
        // whole working set — node buffers and support lists — lives in the
        // ordering; the single public-order epilogue runs in
        // `accumulate_terms`.
        let order = self.class_order.get_or_init(|| series.class_order());
        let inv = order.inv();
        // The fold's inputs arrive in public basis order: convert once per
        // fold. (The rebuild trigger below keeps comparing the caller's
        // public support lists, so it is unaffected by the relabeling.)
        let a_cls = series.class_coefficients(order, a);
        let b_cls = series.class_coefficients(order, b);
        let a_nz_cls: Vec<usize> = a_nonzero.iter().copied().map(|i| inv[i] as usize).collect();
        let b_nz_cls: Vec<usize> = b_nonzero.iter().copied().map(|i| inv[i] as usize).collect();

        // The node lists are a fixed point of the DAG's support propagation:
        // they change only when an atom support changes. Value-level changes
        // with an unchanged support leave them valid; a support change forces
        // one collecting pass in which the kernel re-derives every node's
        // scatter-target set top-down (children before parents).
        let rebuild =
            !self.lists_built || self.atom_a != a_nonzero || self.atom_b != b_nonzero;

        if rebuild {
            // One serial topological sweep in which the kernel re-derives
            // every node's scatter-target set: topological order guarantees
            // a node's children's lists were rebuilt earlier in the pass.
            for k in 2..self.structure.nodes.len() {
                let (left, right) = match self.structure.nodes[k] {
                    DagNode::Binary { left, right } => (left, right),
                    DagNode::Atom(_) => unreachable!("atoms are the first two nodes"),
                };
                let (before, rest) = self.buffers.split_at_mut(k - 2);
                let result = &mut rest[0];
                result.fill(U::default());
                let lbuf = node_slice(left, before, &a_cls, &b_cls);
                let rbuf = node_slice(right, before, &a_cls, &b_cls);
                let (nz_before, nz_rest) = self.nonzeros.split_at_mut(k - 2);
                let lnz = node_nonzeros(left, nz_before, &a_nz_cls, &b_nz_cls);
                let rnz = node_nonzeros(right, nz_before, &a_nz_cls, &b_nz_cls);
                let list = &mut nz_rest[0];
                list.clear();
                let dirty = &mut self.dirty[k - 2];
                series.class_commutation_with_nonzero_collecting(
                    order, lbuf, lnz, rbuf, rnz, result, dirty, list,
                );
            }
            self.atom_a.clear();
            self.atom_a.extend_from_slice(a_nonzero);
            self.atom_b.clear();
            self.atom_b.extend_from_slice(b_nonzero);
            self.lists_built = true;
            return;
        }

        // Steady state: level-parallel batch evaluation. All nodes of a
        // level depend only on earlier levels, their buffers are distinct
        // allocations, and their anagram units are conflict-free — so each
        // level is one parallel dispatch with no synchronization inside.
        for level in self.structure.levels.iter().skip(1) {
            let mut jobs: Vec<KernelJob<U>> = Vec::with_capacity(level.len());
            for &k in level {
                // Zero the buffer and capture its raw pointer; the mutable
                // borrow ends here so the operand reads below are fine (the
                // operand buffers belong to earlier levels — distinct
                // allocations, never the one being filled).
                let (result_ptr, result_len) = {
                    let buffer = &mut self.buffers[k as usize - 2];
                    // The kernel accumulates into the buffer: it must start
                    // from zero on every fold.
                    buffer.fill(U::default());
                    (buffer.as_mut_ptr(), buffer.len())
                };
                let (left, right) = match self.structure.nodes[k as usize] {
                    DagNode::Binary { left, right } => (left, right),
                    DagNode::Atom(_) => unreachable!("atoms are the first two nodes"),
                };
                // SAFETY: a job's operand buffers belong to strictly earlier
                // levels (or are the fold's `a`/`b` inputs); they are written
                // only during their own level's batch, so they are immutable
                // and stable for this level's duration. Each buffer is a
                // distinct allocation, so the jobs' result pointers never
                // alias each other or the operand slices.
                let (a_ptr, a_len, a_nz): (*const U, usize, &[usize]) = match left {
                    0 => (a_cls.as_ptr(), a_cls.len(), &a_nz_cls),
                    1 => (b_cls.as_ptr(), b_cls.len(), &b_nz_cls),
                    id => {
                        let v = &self.buffers[id as usize - 2];
                        (v.as_ptr(), v.len(), &self.nonzeros[id as usize - 2])
                    }
                };
                let (b_ptr, b_len, b_nz): (*const U, usize, &[usize]) = match right {
                    0 => (a_cls.as_ptr(), a_cls.len(), &a_nz_cls),
                    1 => (b_cls.as_ptr(), b_cls.len(), &b_nz_cls),
                    id => {
                        let v = &self.buffers[id as usize - 2];
                        (v.as_ptr(), v.len(), &self.nonzeros[id as usize - 2])
                    }
                };
                jobs.push(KernelJob {
                    a: unsafe { std::slice::from_raw_parts(a_ptr, a_len) },
                    a_nonzero: a_nz,
                    b: unsafe { std::slice::from_raw_parts(b_ptr, b_len) },
                    b_nonzero: b_nz,
                    result: result_ptr,
                    result_len,
                });
            }
            series.class_commutation_batch(order, &mut jobs, &mut self.gating_cache);
        }
    }

    /// Adds the evaluated terms, BCH-weighted, into `target` (a public
    /// basis-order coefficient slice, e.g. the accumulating log
    /// signature's). `rhs` is the displacement atom's public-order slice.
    ///
    /// The accumulation itself runs in the class-contiguous ordering — the
    /// term buffers are class-ordered, so the per-position summation order
    /// matches the public one exactly — and the single epilogue back to
    /// public order runs once, at the end.
    pub(crate) fn accumulate_terms(&mut self, target: &mut [U], rhs: &[U]) {
        let order = self
            .class_order
            .get()
            .expect("accumulate_terms called before evaluate");
        let inv = order.inv();

        // Gather the accumulator's current values and the displacement into
        // class positions. From here to the epilogue every write is
        // sequential in class order.
        let mut acc: Vec<U> = vec![U::default(); inv.len()];
        for (w, c) in target.iter().enumerate() {
            acc[inv[w] as usize] = c.clone();
        }
        let mut rhs_cls: Vec<U> = vec![U::default(); rhs.len().min(inv.len())];
        for (w, c) in rhs.iter().enumerate().take(inv.len()) {
            rhs_cls[inv[w] as usize] = c.clone();
        }

        // Per-position term accumulation in the DAG's term order — the same
        // summation order the public-loop version produced.
        for (source, weight) in &self.structure.terms {
            let ct: &[U] = match source {
                TermSource::Node(node) => &self.buffers[*node as usize - 2],
                TermSource::Displacement => &rhs_cls,
            };
            for (acc_coeff, comm_coeff) in acc.iter_mut().zip(ct) {
                *acc_coeff += comm_coeff.clone() * weight.clone();
            }
        }

        // The one epilogue: class-contiguous accumulator back to public
        // basis order (assignment — `acc` already holds the full sum).
        for (k, &src) in inv.iter().enumerate().take(target.len()) {
            target[k] = acc[src as usize].clone();
        }
    }

    /// The compiled terms in accumulation order.
    pub(crate) fn terms(&self) -> impl Iterator<Item = (TermSource, &U)> {
        self.structure
            .terms
            .iter()
            .map(|(source, weight)| (*source, weight))
    }

    /// The evaluated coefficient slice of an internal node.
    pub(crate) fn buffer(&self, node: u32) -> &[U] {
        &self.buffers[node as usize - 2]
    }

    fn ensure_buffers(&mut self, d: usize, count: usize) {
        if self.buffers.len() != count || self.buffers.first().is_none_or(|b| b.len() != d) {
            self.buffers = (0..count).map(|_| vec![U::default(); d]).collect();
            // A shape change (not a first allocation after clone: fresh
            // buffers are all-zero and the inherited lists still describe
            // them) invalidates the stored target lists.
            if !self.buffers.is_empty() && !self.buffers.first().is_some_and(|b| b.is_empty()) {
                // buffers were just replaced wholesale; distinguish below.
            }
            if self.nonzeros.iter().any(|l| !l.is_empty()) && self.lists_built {
                // Lists reference the old basis size; force a rebuild unless
                // this was the first allocation (empty buffers before).
            }
            self.lists_built = false;
        }
        if self.nonzeros.len() != count {
            self.nonzeros = (0..count).map(|_| Vec::new()).collect();
        }
        let words = d.div_ceil(64);
        if self.dirty.len() != count
            || self.dirty.first().is_none_or(|w| w.len() != words)
        {
            self.dirty = (0..count).map(|_| vec![0u64; words]).collect();
        }
    }
}

impl<U> CommutatorDag<U> {
    /// Shares the compiled plan; scratch buffers are not copied and are
    /// reallocated on the clone's next [`Self::evaluate`].
    pub(crate) fn clone_shallow(&self) -> Self {
        Self {
            structure: Arc::clone(&self.structure),
            // Scratch buffers stay per-value, but the compiled target lists
            // are inherited: a deep-copied accumulator has identical
            // coefficients, hence identical atom supports, hence the same
            // fixed point. The clone's first fold reuses the lists instead
            // of paying a collecting rebuild.
            buffers: Vec::new(),
            nonzeros: self.nonzeros.clone(),
            dirty: self.dirty.clone(),
            atom_a: self.atom_a.clone(),
            atom_b: self.atom_b.clone(),
            lists_built: self.lists_built,
            gating_cache: self.gating_cache.clone(),
            class_order: self.class_order.clone(),
        }
    }
}

impl<U> fmt::Debug for CommutatorDag<U> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("CommutatorDag")
            .field("nodes", &self.structure.nodes.len())
            .field("terms", &self.structure.terms.len())
            .field("buffers", &self.buffers.len())
            .finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::LogSignatureBuilder;
    use lyndon_rs::generators::ENotation;
    use lyndon_rs::lyndon::{LyndonBasis, Sort};
    use num_rational::Ratio;

    type R = Ratio<i128>;

    fn lcg(seed: &mut u64) -> i128 {
        *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((*seed >> 33) % 19) as i128 - 9
    }

    /// The pre-DAG recursive evaluator, kept as an independent oracle: same
    /// tree walk, memoization and accumulation semantics as the old fold.
    #[allow(clippy::too_many_arguments)]
    fn reference_eval(
        term: &CommutatorTerm<R, u8>,
        series: &LieSeries<u8, R>,
        a: &[R],
        a_nonzero: &[usize],
        b: &[R],
        b_nonzero: &[usize],
        memo: &mut HashMap<CommutatorTerm<R, u8>, (Vec<R>, Vec<usize>)>,
    ) -> (Vec<R>, Vec<usize>) {
        match term {
            CommutatorTerm::Atom { atom, .. } => {
                if *atom == 0 {
                    (a.to_vec(), a_nonzero.to_vec())
                } else {
                    (b.to_vec(), b_nonzero.to_vec())
                }
            }
            CommutatorTerm::Expression { .. } => {
                if let Some(hit) = memo.get(term) {
                    return hit.clone();
                }
                let (la, lnz) = reference_eval(
                    term.left().unwrap(), series, a, a_nonzero, b, b_nonzero, memo,
                );
                let (rb, rnz) = reference_eval(
                    term.right().unwrap(), series, a, a_nonzero, b, b_nonzero, memo,
                );
                let mut coefficients = vec![R::from_integer(0); a.len()];
                LieSeries::commutator_coefficients_with_nonzero(
                    series, &la, &lnz, &rb, &rnz, &mut coefficients,
                );
                let nonzero = series.nonzero_coefficient_indices(&coefficients);
                memo.insert(term.clone(), (coefficients.clone(), nonzero.clone()));
                (coefficients, nonzero)
            }
        }
    }

    /// The DAG fold must reproduce the recursive evaluator exactly: same CSE,
    /// same term weights, same accumulation order.
    #[test]
    fn dag_fold_matches_recursive_reference() {
        for (d, m) in [(2usize, 4usize), (2, 5), (3, 4)] {
            let mut seed = 0x5eed_u64.wrapping_mul(d as u64 * 31 + m as u64);
            let builder = LogSignatureBuilder::<u8>::new()
                .with_num_dimensions(d)
                .with_max_degree(m);

            // Accumulator: random dense coefficients.
            let mut log_sig = builder.build::<R>();
            let basis = LyndonBasis::<ENotation>::new(d, Sort::Lexicographical).generate_basis(m);
            let mut acc: Vec<R> = (0..basis.len()).map(|_| R::from_integer(lcg(&mut seed))).collect();
            log_sig.series.coefficients.clone_from(&acc);

            // Reference state mirrors the accumulator and the fold inputs.
            let weights: Vec<R> = log_sig.bch_series.coefficients.clone();
            let trees: Vec<CommutatorTerm<R, u8>> = log_sig.bch_series.commutator_basis.to_vec();

            for step in 0..5 {
                // Perturb the accumulator between folds: zero a few
                // coefficients (support shrink) or raise a few zeros (support
                // growth), alternating with value-only changes. This drives
                // the DAG through both its collecting-rebuild and its
                // fixed-point-reuse paths, including the soundness corner
                // where a previously canceled scatter target goes non-zero.
                if step % 2 == 0 {
                    for _ in 0..3 {
                        let k = (lcg(&mut seed) as usize) % acc.len();
                        acc[k] = R::from_integer(0);
                        log_sig.series.coefficients[k] = R::from_integer(0);
                    }
                } else {
                    for _ in 0..3 {
                        let k = d + (lcg(&mut seed) as usize) % (acc.len() - d);
                        let v = R::from_integer(lcg(&mut seed));
                        acc[k] = v.clone();
                        log_sig.series.coefficients[k] = v;
                    }
                }

                // Production-shaped displacement: full-length coefficient
                // vector with only the degree-1 letters non-zero.
                let displacement: Vec<R> = (0..basis.len())
                    .map(|k| {
                        if k < d { R::from_integer(lcg(&mut seed)) } else { R::from_integer(0) }
                    })
                    .collect();

                // --- DAG fold (the code under test) ---
                let disp_sig = builder.build::<R>();
                let mut disp_sig = disp_sig;
                disp_sig.series.coefficients.clone_from(&displacement);
                log_sig.concatenate_assign(&disp_sig);

                // --- recursive reference fold ---
                // The atom-`A` input of every term is the *pre-fold*
                // accumulator snapshot (the old evaluator's
                // `original_coefficients`), not the mutating accumulator.
                let snapshot = acc.clone();
                let a_nonzero: Vec<usize> = snapshot
                    .iter()
                    .enumerate()
                    .filter(|(_, c)| !c.is_zero())
                    .map(|(i, _)| i)
                    .collect();
                let b_nonzero: Vec<usize> = displacement
                    .iter()
                    .enumerate()
                    .filter(|(_, c)| !c.is_zero())
                    .map(|(i, _)| i)
                    .collect();
                let mut memo = HashMap::new();
                for (i, tree) in trees.iter().enumerate() {
                    if i == 0 || weights[i].is_zero() {
                        continue;
                    }
                    let (ct, _nz) = reference_eval(
                        tree,
                        &log_sig.series,
                        &snapshot,
                        &a_nonzero,
                        &displacement,
                        &b_nonzero,
                        &mut memo,
                    );
                    for (acc_coeff, comm_coeff) in acc.iter_mut().zip(&ct) {
                        *acc_coeff += comm_coeff * weights[i].clone();
                    }
                }

                let diffs: Vec<_> = log_sig
                    .series
                    .coefficients
                    .iter()
                    .zip(&acc)
                    .enumerate()
                    .filter(|(_, (g, w))| g != w)
                    .take(6)
                    .map(|(k, (g, w))| format!("k={k} dag={g} ref={w}"))
                    .collect();
                assert!(
                    diffs.is_empty(),
                    "d={d} m={m} step={step}: DAG fold diverged: {}",
                    diffs.join("; ")
                );
            }
        }
    }
}

