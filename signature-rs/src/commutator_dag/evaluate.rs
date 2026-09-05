//! Single-fold plan evaluation, batch eligibility and the public-order
//! term-accumulation epilogue.

use lie_rs::{ClassOrderedCommutation, KernelJob, LieSeries};
use lyndon_rs::generators::Generator;
use num_traits::{One, Zero};
use std::hash::Hash;
use std::ops::{AddAssign, MulAssign, Neg};
use std::sync::Arc;

use super::structure::{DagNode, TermSource};
use super::{CommutatorDag, node_nonzeros, node_slice};

impl<U> CommutatorDag<U>
where
    U: Clone
        + Default
        + One
        + Zero
        + Eq
        + MulAssign
        + Neg<Output = U>
        + Hash
        + AddAssign
        + 'static,
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
        let rebuild = !self.lists_built || self.atom_a != a_nonzero || self.atom_b != b_nonzero;

        if rebuild {
            // Copy-on-write: a dag sharing an adopted snapshot's lists must
            // never mutate them through the `Arc`. Unwrap when exclusively
            // held (free), clone when shared, mutate the owned copy, and
            // install it as this dag's own `Arc`.
            let mut nz = Arc::unwrap_or_clone(std::mem::take(&mut self.nonzeros));
            let mut dt = Arc::unwrap_or_clone(std::mem::take(&mut self.dirty));
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
                let (nz_before, nz_rest) = nz.split_at_mut(k - 2);
                let lnz = node_nonzeros(left, nz_before, &a_nz_cls, &b_nz_cls);
                let rnz = node_nonzeros(right, nz_before, &a_nz_cls, &b_nz_cls);
                let list = &mut nz_rest[0];
                list.clear();
                let dirty = &mut dt[k - 2];
                series.class_commutation_with_nonzero_collecting(
                    order, lbuf, lnz, rbuf, rnz, result, dirty, list,
                );
            }
            self.nonzeros = Arc::new(nz);
            self.dirty = Arc::new(dt);
            self.atom_a.clear();
            self.atom_a.extend_from_slice(a_nonzero);
            self.atom_b.clear();
            self.atom_b.extend_from_slice(b_nonzero);
            self.lists_built = true;
            return;
        }

        // Steady state: one fold-parallel dispatch. All levels' jobs are
        // built up front (zeroing every buffer first — the kernel
        // accumulates); the engine plans class-contiguous balanced packs per
        // level and walks them with dynamic pack claims, the level counters
        // ordering cross-level operand reads.
        let mut levels: Vec<Vec<KernelJob<U>>> = self
            .structure
            .levels
            .iter()
            .skip(1)
            .map(|level| {
                level
                    .iter()
                    .filter_map(|&k| match self.structure.nodes[k as usize] {
                        DagNode::Binary { .. } => {
                            // The kernel accumulates into the buffer: it must
                            // start from zero on every fold.
                            let buffer = &mut self.buffers[k as usize - 2];
                            buffer.fill(U::default());
                            Some(KernelJob {
                                a: std::ptr::null(),
                                a_len: 0,
                                a_nonzero: &[],
                                b: std::ptr::null(),
                                b_len: 0,
                                b_nonzero: &[],
                                result: buffer.as_mut_ptr(),
                                result_len: buffer.len(),
                                a_shift: lie_rs::IDENTITY_SHIFTS.as_ptr(),
                                b_shift: lie_rs::IDENTITY_SHIFTS.as_ptr(),
                                r_shift: lie_rs::IDENTITY_SHIFTS.as_ptr(),
                            })
                        }
                        DagNode::Atom(_) => None,
                    })
                    .collect()
            })
            .collect();
        // Fill the operand views after the zeroing pass (the mutable borrows
        // above have ended; operand buffers belong to strictly earlier
        // levels or are the fold's converted inputs, and each result buffer
        // is a distinct allocation, so the raw views never alias within a
        // level).
        for (li, level) in self.structure.levels.iter().enumerate().skip(1) {
            for (jj, &k) in level.iter().enumerate() {
                let (left, right) = match self.structure.nodes[k as usize] {
                    DagNode::Binary { left, right } => (left, right),
                    DagNode::Atom(_) => continue,
                };
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
                // SAFETY: the operand buffers are stable and unmodified for
                // this level's duration (they belong to earlier levels,
                // ordered by the engine's level counters); the result
                // pointers of one level never alias each other.
                let job = &mut levels[li - 1][jj];
                job.a = a_ptr;
                job.a_len = a_len;
                job.a_nonzero = a_nz;
                job.b = b_ptr;
                job.b_len = b_len;
                job.b_nonzero = b_nz;
            }
        }
        series.class_commutation_batch_fold(
            order,
            &mut levels,
            &mut self.gating_cache,
        );
    }

    /// Batch-readiness: the node lists are the built fixed point for these
    /// atom supports, and the accumulator's support already equals the
    /// reachable set — the union of every node's scatter-target list and
    /// the displacement support. Future accumulation can then only change
    /// values at already-nonzero positions (a later cancellation merely
    /// shrinks the support, and the superset lists stay sound), so a whole
    /// run of folds can share one plan and one set of gating masks.
    pub(crate) fn batch_eligible(&self, a_nonzero: &[usize], b_nonzero: &[usize]) -> bool {
        if !self.lists_built || self.atom_a != a_nonzero || self.atom_b != b_nonzero {
            return false;
        }
        let Some(order) = self.class_order.get() else {
            return false;
        };
        let inv = order.inv();
        let words = inv.len();
        let words64 = words.div_ceil(64);
        let mut reach = vec![0u64; words64];
        for list in self.nonzeros.iter() {
            for &i in list {
                reach[i / 64] |= 1 << (i % 64);
            }
        }
        for &w in b_nonzero {
            let i = inv[w] as usize;
            reach[i / 64] |= 1 << (i % 64);
        }
        let mut acc = vec![0u64; words64];
        for &w in a_nonzero {
            let i = inv[w] as usize;
            acc[i / 64] |= 1 << (i % 64);
        }
        acc == reach
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
        // summation order the public-loop version produced. `raw_mul`/
        // `raw_add_assign` skip the float wrappers' per-op NaN checks
        // (raw-float fast path, see lie-rs `raw_mul`'s NaN policy).
        for (source, weight) in self.structure.terms.iter() {
            let ct: &[U] = match source {
                TermSource::Node(node) => &self.buffers[*node as usize - 2],
                TermSource::Displacement => &rhs_cls,
            };
            for (acc_coeff, comm_coeff) in acc.iter_mut().zip(ct) {
                lie_rs::raw_add_assign(acc_coeff, &lie_rs::raw_mul(comm_coeff, weight));
            }
        }

        // The one epilogue: class-contiguous accumulator back to public
        // basis order (assignment — `acc` already holds the full sum).
        for (k, &src) in inv.iter().enumerate().take(target.len()) {
            target[k] = acc[src as usize].clone();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use commutator_rs::CommutatorTerm;
    use crate::LogSignatureBuilder;
    use crate::commutator_dag::DagStructure;
    use lie_rs::GatingCache;
    use lyndon_rs::generators::ENotation;
    use lyndon_rs::lyndon::{LyndonBasis, Sort};
    use num_rational::Ratio;
    use num_traits::{One, Zero};
    use std::collections::HashMap;
    use std::sync::OnceLock;

    type R = Ratio<i128>;

    fn lcg(seed: &mut u64) -> i128 {
        *seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*seed >> 33) % 19) as i128 - 9
    }

    fn atom(atom: u8) -> CommutatorTerm<R, u8> {
        CommutatorTerm::Atom {
            coefficient: R::one(),
            atom,
        }
    }

    fn bracket(left: CommutatorTerm<R, u8>, right: CommutatorTerm<R, u8>) -> CommutatorTerm<R, u8> {
        let degree = left.degree() + 1;
        CommutatorTerm::Expression {
            coefficient: R::one(),
            left: Box::new(left),
            right: Box::new(right),
            degree,
        }
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
                    term.left().unwrap(),
                    series,
                    a,
                    a_nonzero,
                    b,
                    b_nonzero,
                    memo,
                );
                let (rb, rnz) = reference_eval(
                    term.right().unwrap(),
                    series,
                    a,
                    a_nonzero,
                    b,
                    b_nonzero,
                    memo,
                );
                let mut coefficients = vec![R::from_integer(0); a.len()];
                LieSeries::commutator_coefficients_with_nonzero(
                    series,
                    &la,
                    &lnz,
                    &rb,
                    &rnz,
                    &mut coefficients,
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
            let mut acc: Vec<R> = (0..basis.len())
                .map(|_| R::from_integer(lcg(&mut seed)))
                .collect();
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
                        if k < d {
                            R::from_integer(lcg(&mut seed))
                        } else {
                            R::from_integer(0)
                        }
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

    /// Liveness of the one-dispatch fold under nested parallelism. When the
    /// shared pool is saturated by outer tasks, a fold's slot tasks queue
    /// behind them; a fold that blocked on a *queued* releaser's pack would
    /// never finish (the fixed slot→pack assignment starved exactly there).
    /// Pack claims go to whoever runs, so folds nested inside an
    /// oversubscribed outer dispatch must complete — and match the serial
    /// result bit for bit, since a claim relabels which slot sweeps a pack,
    /// never the order of adds within a word.
    #[test]
    fn dag_fold_survives_nested_parallel_oversubscription() {
        use rayon::iter::{IntoParallelIterator, ParallelIterator};
        use std::sync::atomic::AtomicBool;

        let d = 3usize;
        let m = 6usize;
        let mut seed = 0x51ea_u64;
        let builder = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(d)
            .with_max_degree(m);
        let basis = LyndonBasis::<ENotation>::new(d, Sort::Lexicographical).generate_basis(m);

        // One accumulator and one displacement, folded once on the calling
        // thread to produce the expected post-fold coefficients.
        let acc: Vec<R> = (0..basis.len())
            .map(|_| R::from_integer(lcg(&mut seed)))
            .collect();
        let displacement: Vec<R> = (0..basis.len())
            .map(|k| {
                if k < d {
                    R::from_integer(lcg(&mut seed))
                } else {
                    R::from_integer(0)
                }
            })
            .collect();
        let mut template = builder.build::<R>();
        template.series.coefficients.clone_from(&acc);
        let mut disp_sig = builder.build::<R>();
        disp_sig.series.coefficients.clone_from(&displacement);
        let mut reference = template.clone();
        reference.concatenate_assign(&disp_sig);
        let expected = reference.series.coefficients.clone();

        // Watchdog: a liveness bug shows up as a hang, not a wrong value —
        // turn it into a hard failure (abort fails the test process).
        let done = std::sync::Arc::new(AtomicBool::new(false));
        let flag = done.clone();
        std::thread::spawn(move || {
            for _ in 0..600 {
                std::thread::sleep(std::time::Duration::from_millis(100));
                if flag.load(std::sync::atomic::Ordering::SeqCst) {
                    return;
                }
            }
            eprintln!("nested fold stress did not complete in 60s: fold walk starved");
            std::process::abort();
        });

        // 2× oversubscription: 64 nested folds on the default pool.
        (0..64usize).into_par_iter().for_each(|_| {
            let mut sig = template.clone();
            sig.concatenate_assign(&disp_sig);
            assert_eq!(sig.series.coefficients, expected);
        });
        done.store(true, std::sync::atomic::Ordering::SeqCst);
    }

    /// The stripped hand DAG is executable: evaluating the compacted
    /// 4-node plan (atoms + [A,B] + [[A,B],A]) through the production
    /// rebuild path must match the recursive oracle bit for bit.
    #[test]
    fn strip_then_evaluate_matches_recursive_oracle() {
        let d = 3usize;
        let m = 3usize;
        let mut seed = 0xdead_u64;
        let builder = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(d)
            .with_max_degree(m);
        let basis = LyndonBasis::<ENotation>::new(d, Sort::Lexicographical).generate_basis(m);

        // The compacted DAG from the renumbering test: [A,B] and [[A,B],A].
        let mut dag = CommutatorDag::<R> {
            structure: std::sync::Arc::new(DagStructure {
                nodes: vec![
                    DagNode::Atom(0),
                    DagNode::Atom(1),
                    DagNode::Binary {
                        left: 0,
                        right: 1,
                    },
                    DagNode::Binary {
                        left: 2,
                        right: 0,
                    },
                ],
                levels: vec![vec![0, 1], vec![2], vec![3]],
                terms: vec![
                    (TermSource::Node(2), R::from_integer(3)),
                    (TermSource::Node(3), R::from_integer(-7)),
                ],
            }),
            buffers: Vec::new(),
            nonzeros: Arc::new(Vec::new()),
            dirty: Arc::new(Vec::new()),
            atom_a: Vec::new(),
            atom_b: Vec::new(),
            lists_built: false,
            gating_cache: GatingCache::default(),
            class_order: OnceLock::new(),
        };

        let a: Vec<R> = (0..basis.len())
            .map(|_| R::from_integer(lcg(&mut seed)))
            .collect();
        let b: Vec<R> = (0..basis.len())
            .map(|k| {
                if k < d {
                    R::from_integer(lcg(&mut seed))
                } else {
                    R::from_integer(0)
                }
            })
            .collect();
        let a_nonzero: Vec<usize> = a
            .iter()
            .enumerate()
            .filter(|(_, c)| !c.is_zero())
            .map(|(i, _)| i)
            .collect();
        let b_nonzero: Vec<usize> = b
            .iter()
            .enumerate()
            .filter(|(_, c)| !c.is_zero())
            .map(|(i, _)| i)
            .collect();

        dag.evaluate(&builder.build::<R>().series, &a, &a_nonzero, &b, &b_nonzero);

        // Recursive oracle on the equivalent bracket trees. Node buffers
        // live in the DAG's class-contiguous order (the public-order
        // epilogue only runs in `accumulate_terms`), so the oracle's
        // public-order coefficients are gathered through the same
        // `class[inv[w]] = public[w]` map the fold's inputs use before
        // comparing.
        let t2 = bracket(atom(0), atom(1)); // node 2 = [A, B]
        let t3 = bracket(t2.clone(), atom(0)); // node 3 = [[A,B], A]
        let series = builder.build::<R>().series;
        let inv = dag
            .class_order
            .get()
            .expect("evaluate built the class order")
            .inv();
        let mut memo = HashMap::new();
        for (node, tree) in [(2u32, &t2), (3, &t3)] {
            let (ct, _nz) = reference_eval(
                tree,
                &series,
                &a,
                &a_nonzero,
                &b,
                &b_nonzero,
                &mut memo,
            );
            let mut cls = vec![R::zero(); ct.len()];
            for (w, value) in ct.iter().enumerate() {
                cls[inv[w] as usize] = value.clone();
            }
            assert_eq!(
                dag.buffers[node as usize - 2], cls,
                "node {node} buffer diverged from the recursive oracle"
            );
        }
    }
}
