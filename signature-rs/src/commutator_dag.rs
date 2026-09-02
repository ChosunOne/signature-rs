//! Compiled evaluation plan for the BCH concatenation expression: shared
//! subtrees interned to node ids, evaluated as a straight-line loop over
//! reused per-node buffers instead of recursive per-fold evaluation.

use std::cell::UnsafeCell;
use std::collections::HashMap;
use std::fmt;
use std::hash::Hash;
use std::ops::{AddAssign, MulAssign, Neg};
use std::sync::{Arc, OnceLock};

use commutator_rs::CommutatorTerm;
use lie_rs::{
    ClassBatchStage, ClassOrder, ClassOrderedCommutation, GatingCache, KernelJob, LieSeries,
    plan_class_sweep_stages, planned_sweep_entries, run_class_batch, work_adaptive_slots,
};
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

/// Raw-pointer wrapper for the batch's cross-stage working buffers. All
/// access goes through these pointers (never through a `&mut` taken after
/// derivation); disjointness comes from the block ranges within a stage,
/// and the batch's stage counters order every cross-stage dependency.
struct SendRaw<T>(*mut T);
unsafe impl<T> Send for SendRaw<T> {}
unsafe impl<T> Sync for SendRaw<T> {}

/// Raw-pointer wrapper for read-only cross-stage inputs (the per-fold
/// displacement slices).
struct SendConst<T>(*const T);
unsafe impl<T> Send for SendConst<T> {}
unsafe impl<T> Sync for SendConst<T> {}

/// The node buffers' `(pointer, length)` table for the batch's block
/// closures: the pointers are raw, so the table itself needs the Send/Sync
/// promise (same argument as [`SendRaw`]).
struct SendBufPtrs<U>(Vec<(*mut U, usize)>);
unsafe impl<U> Send for SendBufPtrs<U> {}
unsafe impl<U> Sync for SendBufPtrs<U> {}

/// One accumulate run of a block's work list: a contiguous (class
/// position, source slot) stretch of one BCH term whose positions fall in
/// the block's range.
struct AccumRun<U> {
    /// Index into `DagStructure::terms` for the run's BCH weight.
    term: u32,
    /// First class position of the run.
    g0: u32,
    /// Number of class positions in the run.
    len: u32,
    /// Source pointer: the node's compact slot holding position `g0` (the
    /// fold displacement at `g0` for the displacement term); advances with
    /// `g0`. Mutable because fused runs write `U::default()` back after
    /// consuming each slot (the pointer is derived from the batch's
    /// mutable compact/displacement buffers).
    src: *mut U,
    /// Fused zero-after-add: each source slot is written back to
    /// `U::default()` as it is consumed. Set only for a node referenced by
    /// exactly one term (its slots are dead once that term has read them);
    /// the displacement term's `b_cls` source is rewritten wholesale by
    /// the next fold's gather, so it never needs this.
    zero: bool,
}

/// One zero run of a block's work list: a contiguous stretch of one
/// compact node buffer's slots inside the block's class range.
struct ZeroRun<U> {
    /// Target compact buffer pointer.
    ptr: *mut U,
    /// First compact slot of the run.
    s0: u32,
    /// Number of compact slots in the run.
    len: u32,
}

/// A block's per-batch work lists: the accumulate runs and zero runs
/// intersecting the block's class-position range (see `fold_batch`).
struct BlockWork<U> {
    accum: Vec<AccumRun<U>>,
    zero: Vec<ZeroRun<U>>,
}

/// The per-block work lists of one batch, shared by the batch's block
/// closures through an `Arc` (raw-pointer bearing, same Send/Sync argument
/// as [`SendBufPtrs`]: the pointers target stable frame-owned allocations
/// whose cross-stage ordering the stage counters provide).
struct SendBlockWork<U>(std::sync::Arc<Vec<BlockWork<U>>>);
unsafe impl<U> Send for SendBlockWork<U> {}
unsafe impl<U> Sync for SendBlockWork<U> {}

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

/// Groups a sorted ascending position list into maximal runs of consecutive
/// positions within one block (position `p` belongs to block
/// `p * blocks / d`), invoking `emit(block, first_position, run_length)`
/// per run in ascending order. Consecutive positions of one block merge
/// into a single run; a block boundary or a position gap starts a new one.
fn for_each_position_run(
    positions: &[u32],
    blocks: usize,
    d: usize,
    mut emit: impl FnMut(usize, u32, u32),
) {
    let mut i = 0;
    while i < positions.len() {
        let first = positions[i];
        let block = first as usize * blocks / d;
        let mut last = first;
        let mut j = i + 1;
        while j < positions.len()
            && positions[j] == last + 1
            && positions[j] as usize * blocks / d == block
        {
            last = positions[j];
            j += 1;
        }
        emit(block, first, last - first + 1);
        i = j;
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
                    let left = intern(
                        term.left().expect("expression has a left child"),
                        nodes,
                        interned,
                    );
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
            structure: Arc::new(DagStructure {
                nodes,
                terms,
                levels,
            }),
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
        series.class_commutation_batch_fold(order, &mut levels, &mut self.gating_cache);
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
        for list in &self.nonzeros {
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

    /// Folds `rhss` (public-order displacement slices) into the class-space
    /// accumulator `acc` (public order, in-out) as ONE continuous stage
    /// chain: per fold, a gather stage brings the displacement into class
    /// space, the DAG's levels sweep (reading the class-space accumulator
    /// and the shared displacement buffer), and an accumulate stage adds
    /// every term into the class-space accumulator, zeroing the node
    /// buffers afterwards for the next fold. One plan and one job table
    /// serve the whole batch; the accumulator never leaves class space, so
    /// the per-fold public→class→public round-trip disappears and the
    /// plan/gating cost is paid once, not per fold.
    ///
    /// Contract (checked by the driver): the accumulator's support equals
    /// the reachable set (the union of the node lists and the displacement
    /// support, see [`Self::batch_eligible`]) and all displacements share
    /// the support `b_nonzero`. Level-0 gating then uses masks that cover
    /// every position any later fold can touch, which stays sound even if
    /// mid-batch accumulation cancels values to zero — the node lists are
    /// scatter-target supersets, so a shrinking support only loses
    /// tightness, never targets.
    pub(crate) fn fold_batch<T>(
        &mut self,
        series: &mut LieSeries<T, U>,
        rhss: &[&[U]],
        a_nonzero: &[usize],
        b_nonzero: &[usize],
    ) where
        T: Clone + Ord + Generator + Hash + Eq,
        U: Send + Sync,
    {
        let internal = self.structure.nodes.len() - 2;
        self.ensure_buffers(series.coefficients.len(), internal);
        let order = self
            .class_order
            .get()
            .expect("batch fold requires the class order of a prior steady fold");
        let inv = order.inv();
        let perm = order.perm();
        let d = series.coefficients.len();

        // The batch's class-space working state. The accumulator stays
        // class-ordered across folds: fold f's sweep reads it as its `a`
        // operand and fold f's accumulate adds the terms in place. Both
        // buffers are written by block stages and read by the sweeps'
        // operand slices; all access goes through the cells' raw pointers
        // (the UnsafeCell makes that interleaving defined) and the stage
        // counters order every write-before-read across stages.
        let mut acc_cls = vec![U::default(); d];
        for (w, c) in series.coefficients.iter().enumerate().take(d) {
            acc_cls[inv[w] as usize] = c.clone();
        }
        let acc_cell = UnsafeCell::new(acc_cls);
        let b_cell = UnsafeCell::new(vec![U::default(); d]);
        // SAFETY: the cells are not moved for the dispatch's duration, so
        // the heap buffers behind their Vecs are stable; everything below
        // goes through these data pointers (the sweeps read shared slices
        // derived from them, the block stages read and write through
        // them), and the stage counters order every write-before-read.
        let acc_data: *mut U = unsafe { (*acc_cell.get()).as_mut_ptr() };
        let b_data: *mut U = unsafe { (*b_cell.get()).as_mut_ptr() };

        // Compact per-node scratch buffers. A node's sweep only ever
        // scatters into the degree slices its recorded support touches
        // (averaging ~1.5 slices per node), and the class space is
        // degree-grouped, so storing whole slices costs ~4-6x less memory
        // than one full-d buffer per node — deep DAGs (hundreds of nodes)
        // drop from megabytes to tens/hundreds of KB and the sweeps' scatter
        // RMWs become L2-resident. Layout per node: its active degree
        // slices concatenated in degree order; class position `x` of degree
        // `δ` lives at `x - shift[δ]` (the per-degree shift table the jobs
        // carry; a full-d buffer uses the all-zero identity table).
        // Per-class-position Lyndon degrees (the class space is
        // degree-grouped, so this is non-decreasing). NOTE: class positions,
        // not the public basis order — the batch works entirely in the
        // class-contiguous space.
        let pos_degree: Vec<u8> = order.degree_cls().to_vec();
        let max_deg = pos_degree.iter().copied().max().unwrap_or(0) as usize;
        // Degree-slice starts: the class space is degree-grouped, so slice
        // `dd` = [first(dd), first(dd+1)) where first(x) is the first class
        // position of degree x (d when degree x is empty; backfilled to the
        // next populated degree's start so empty slices are zero-length).
        let mut deg_start: Vec<u32> = vec![d as u32; max_deg + 2];
        for (pos, &dg) in pos_degree.iter().enumerate() {
            let dg = dg as usize;
            if (deg_start[dg] as usize) > pos {
                deg_start[dg] = pos as u32;
            }
        }
        deg_start[max_deg + 1] = d as u32;
        for dd in (0..=max_deg).rev() {
            if deg_start[dd] > deg_start[dd + 1] {
                deg_start[dd] = deg_start[dd + 1];
            }
        }

        // Identity shifts for the full-d atom buffers (acc_cls / b_cls).
        // The compact per-node buffers are built after the sweep planner
        // returns the jobs' exact batch scatter sets (the recorded node
        // lists are only a union-level fixed point — the batch eligibility
        // checks the union, not each list — so a level-0 job's list,
        // recorded under an earlier, smaller accumulator support, can miss
        // positions the batch's first fold scatters; sizing from the
        // recorded lists would be unsound).
        let zero_shifts = vec![0u32; max_deg + 2];
        let zero_shifts_ptr: *const u32 = zero_shifts.as_ptr();

        // One job table for the whole batch: level-0 jobs read the
        // class-space accumulator and displacement buffers (fixed
        // allocations, values refreshed in-graph by the gather stages)
        // with the reachable/common supports; interior jobs read child
        // node buffers with the children's scatter-target lists — exactly
        // the per-fold operand fill. The jobs' `result`/`r_shift` are
        // wired after the planner returns the exact per-node scatter sets.
        let a_nz_cls: Vec<usize> = a_nonzero.iter().map(|&i| inv[i] as usize).collect();
        let b_nz_cls: Vec<usize> = b_nonzero.iter().map(|&i| inv[i] as usize).collect();
        let mut levels_jobs: Vec<Vec<KernelJob<U>>> = self
            .structure
            .levels
            .iter()
            .skip(1)
            .map(|level| {
                level
                    .iter()
                    .filter_map(|&k| match self.structure.nodes[k as usize] {
                        DagNode::Binary { .. } => Some(KernelJob {
                            a: acc_data,
                            a_len: d,
                            a_nonzero: &a_nz_cls,
                            b: b_data,
                            b_len: d,
                            b_nonzero: &b_nz_cls,
                            result: std::ptr::null_mut(),
                            result_len: 0,
                            a_shift: zero_shifts_ptr,
                            b_shift: zero_shifts_ptr,
                            r_shift: zero_shifts_ptr,
                        }),
                        DagNode::Atom(_) => None,
                    })
                    .collect()
            })
            .collect();
        // SAFETY: the operand slices are derived from the cell data
        // pointers and the node buffers (stable allocations, no resize
        // during the dispatch) and only read (by the sweeps); the block
        // stages write through the same raws between stages — the
        // UnsafeCell makes that interleaving defined, and the stage
        // counters order it. The allocations live to the end of this
        // function (the dispatch joins before they drop).
        // Local copies of the node lists for the jobs' operand views (the
        // lists are re-recorded from the planner's exact scatter sets
        // below, so the jobs must not borrow `self.nonzeros`).
        let lists_local: Vec<Vec<usize>> = self.nonzeros.clone();
        for (li, level) in self.structure.levels.iter().enumerate().skip(1) {
            for (jj, &k) in level.iter().enumerate() {
                let (left, right) = match self.structure.nodes[k as usize] {
                    DagNode::Binary { left, right } => (left, right),
                    DagNode::Atom(_) => continue,
                };
                let (a_ptr, a_len, a_sh, a_nz): (*const U, usize, *const u32, &[usize]) = match left
                {
                    0 => (acc_data as *const U, d, zero_shifts_ptr, &a_nz_cls),
                    1 => (b_data as *const U, d, zero_shifts_ptr, &b_nz_cls),
                    id => (
                        std::ptr::null(),
                        0,
                        zero_shifts_ptr,
                        &lists_local[id as usize - 2],
                    ),
                };
                let (b_ptr, b_len, b_sh, b_nz): (*const U, usize, *const u32, &[usize]) =
                    match right {
                        0 => (acc_data as *const U, d, zero_shifts_ptr, &a_nz_cls),
                        1 => (b_data as *const U, d, zero_shifts_ptr, &b_nz_cls),
                        id => (
                            std::ptr::null(),
                            0,
                            zero_shifts_ptr,
                            &lists_local[id as usize - 2],
                        ),
                    };
                let job = &mut levels_jobs[li - 1][jj];
                job.a = a_ptr;
                job.a_len = a_len;
                job.a_nonzero = a_nz;
                job.b = b_ptr;
                job.b_len = b_len;
                job.b_nonzero = b_nz;
                job.a_shift = a_sh;
                job.b_shift = b_sh;
                job.r_shift = zero_shifts_ptr;
            }
        }

        // Plan once to get the jobs' exact batch scatter sets (the gating
        // depends only on the fixed support lists, so this plan's stages
        // are throwaway — they borrow the placeholder jobs and are dropped
        // before the jobs are rewired; the final plan below hits the
        // gating cache).
        let (_, scatter_sets) =
            plan_class_sweep_stages(series, order, &levels_jobs, &mut self.gating_cache, true);

        // Record the true scatter sets as the node lists: they bound this
        // batch's sweeps exactly, keep the union-level eligibility fixed
        // point for later batches, and stay sound supersets for the per-
        // fold path's gating.
        {
            let mut updated = std::mem::take(&mut self.nonzeros);
            for (li, level) in self.structure.levels.iter().enumerate().skip(1) {
                for (jj, &k) in level.iter().enumerate() {
                    if matches!(self.structure.nodes[k as usize], DagNode::Binary { .. }) {
                        updated[k as usize - 2] = scatter_sets[li - 1][jj]
                            .iter()
                            .map(|&p| p as usize)
                            .collect();
                    }
                }
            }
            self.nonzeros = updated;
        }

        // Compact per-node scratch buffers, sized from the exact scatter
        // sets: per node, its support's degree slices concatenated in
        // degree order; class position `x` of degree `δ` lives at
        // `x - shift[δ]`.
        struct CompactSlice {
            base: u32,
            rs: u32,
            re: u32,
        }
        // Indexed BY NODE ID (`k - 2`), not by (level, job) position — the
        // DAG's levels are not id-ordered (a late-interned shallow node can
        // share a level with early deep nodes), so (level, job) iteration
        // order must not be conflated with node identity.
        let internal = self.buffers.len();
        let mut compact: Vec<Vec<U>> = Vec::with_capacity(internal);
        let mut shifts: Vec<Vec<u32>> = Vec::with_capacity(internal);
        let mut layouts: Vec<Vec<CompactSlice>> = Vec::with_capacity(internal);
        for _ in 0..internal {
            compact.push(Vec::new());
            shifts.push(vec![0u32; max_deg + 2]);
            layouts.push(Vec::new());
        }
        for (li, level) in self.structure.levels.iter().enumerate().skip(1) {
            for (jj, &k) in level.iter().enumerate() {
                let set: &[u32] =
                    if matches!(self.structure.nodes[k as usize], DagNode::Binary { .. }) {
                        &scatter_sets[li - 1][jj]
                    } else {
                        &[]
                    };
                let ki = k as usize - 2;
                let mut degs: Vec<u8> = set.iter().map(|&p| pos_degree[p as usize]).collect();
                degs.sort_unstable();
                degs.dedup();
                let mut buf: Vec<U> = Vec::new();
                let shift = &mut shifts[ki];
                let slices = &mut layouts[ki];
                for dg in degs {
                    let rs = deg_start[dg as usize];
                    let re = deg_start[dg as usize + 1];
                    let base = buf.len() as u32;
                    buf.resize(buf.len() + (re - rs) as usize, U::default());
                    shift[dg as usize] = rs - base;
                    slices.push(CompactSlice { base, rs, re });
                }
                compact[ki] = buf;
            }
        }
        let compact_ptrs: Vec<(*mut U, usize)> = compact
            .iter_mut()
            .map(|b| (b.as_mut_ptr(), b.len()))
            .collect();

        // Wire the jobs' results and the interior jobs' operand views to
        // the compact buffers (the nodes' slices cover their recorded
        // supports, which is exactly what the presence tests load).
        for (li, level) in self.structure.levels.iter().enumerate().skip(1) {
            for (jj, &k) in level.iter().enumerate() {
                let job = &mut levels_jobs[li - 1][jj];
                let (rp, rl) = compact_ptrs[k as usize - 2];
                job.result = rp;
                job.result_len = rl;
                job.r_shift = shifts[k as usize - 2].as_ptr();
                let (left, right) = match self.structure.nodes[k as usize] {
                    DagNode::Binary { left, right } => (left, right),
                    DagNode::Atom(_) => continue,
                };
                if left > 1 {
                    let (p, l) = compact_ptrs[left as usize - 2];
                    job.a = p as *const U;
                    job.a_len = l;
                    job.a_shift = shifts[left as usize - 2].as_ptr();
                }
                if right > 1 {
                    let (p, l) = compact_ptrs[right as usize - 2];
                    job.b = p as *const U;
                    job.b_len = l;
                    job.b_shift = shifts[right as usize - 2].as_ptr();
                }
            }
        }

        // Final plan with the wired jobs: the gating cache makes this cheap
        // (the support lists are unchanged); the stages reference the final
        // job table the sweeps read through.
        let (sweep_stages, _) =
            plan_class_sweep_stages(series, order, &levels_jobs, &mut self.gating_cache, true);

        // Block ranges over class positions. The block count derives from
        // the SAME work-adaptive policy the walk will use (below): blocks
        // exist to spread the gather/accumulate phases across slots, so
        // cutting more blocks than the policy's slot count only multiplies
        // the walk's per-pack claim/publish protocol (78 micro-blocks per
        // stage × ~1000 folds of atomic RMWs dwarfed the gather itself and
        // kept the tiny-grid e2e 2.3x above its 1t floor even at slots=1).
        // The policy sees the per-fold sweep work; the two block stages a
        // fold adds are counted at their post-cut element count, which for
        // the slot decision changes nothing (both terms are « QUANTUM
        // exactly when the sweep term is).
        let threads = rayon::current_num_threads().max(1);
        let sweep_refs: Vec<&ClassBatchStage<U>> = sweep_stages.iter().collect();
        let slots = work_adaptive_slots(threads, d.max(1), planned_sweep_entries(&sweep_refs));
        let blocks = d.min(4 * slots).max(1);
        let ranges: Vec<(u32, u32)> = (0..blocks)
            .map(|b| ((d * b / blocks) as u32, (d * (b + 1) / blocks) as u32))
            .collect();

        // Block closures, kept alive for the whole dispatch. Z zeroes the
        // compact node buffers before the first in-batch sweep (the
        // accumulate stages then maintain the zeroing); per fold, BG gathers
        // the displacement into class space and C accumulates the terms.
        //
        // Those runs are materialized once per batch, here at plan time —
        // and only over each source's LIVE slots. A node's sweep scatters
        // exactly onto its recorded scatter set (the planner's exact write
        // set, re-recorded into `self.nonzeros` above), so compact slots
        // outside it are padding that no sweep ever touches — permanently
        // `U::default()`. Their `weight × (+0.0)` adds are value-zero and
        // dropped: the only observable difference is the sign bit of an
        // exact-zero accumulator slot (adding +0.0 turns a −0.0 into +0.0;
        // adding −0.0 is a bitwise identity), never the `==` value. The
        // runs keep the full walk's per-position add order (ascending term
        // index; within a term, ascending degree slice and ascending
        // position), so every position's add sequence — value, weight, term
        // order — is unchanged. Zeroing is fused into the single
        // referencing term's runs (a node is the canonical root of at most
        // one BCH term — distinct Lyndon words have distinct canonical
        // bracketings — so no other run reads the slot after that term's
        // add); nodes no term accumulates (shared interior brackets, or
        // defensively a multi-rooted node) keep explicit zero runs over
        // their live slots.
        let mut block_work: Vec<BlockWork<U>> = (0..blocks)
            .map(|_| BlockWork {
                accum: Vec::new(),
                zero: Vec::new(),
            })
            .collect();
        // Referencing terms per node: count 1 → that term's accum runs fuse
        // the slot zeroing; count 0 (interior node) or ≥ 2 (defensive) →
        // explicit zero runs.
        let mut node_terms: Vec<Vec<u32>> = vec![Vec::new(); internal];
        for (ti, (source, _)) in self.structure.terms.iter().enumerate() {
            if let TermSource::Node(k) = source {
                node_terms[*k as usize - 2].push(ti as u32);
            }
        }
        for (ti, (source, _)) in self.structure.terms.iter().enumerate() {
            match source {
                TermSource::Displacement => {
                    // The gather rewrites every b_cls slot each fold with
                    // `rhs[perm[g]]` (`U::default()` beyond the rhs
                    // length), and the batch driver guarantees every rhs is
                    // value-zero outside the common support over the
                    // kernel-reachable positions
                    // `[0, degree_start(max_degree))` — exactly the prefix
                    // below `tail` (class and public layouts share degree-
                    // slice boundaries). The degree-`max_degree` tail lies
                    // beyond that check, so its whole contiguous range
                    // stays covered: its gather-written values fold
                    // unconditionally, exactly as before. Everything
                    // skipped in the prefix is a value-zero add (see the
                    // sign-of-zero note above).
                    let tail = deg_start[max_deg];
                    let mut sup: Vec<u32> = b_nz_cls
                        .iter()
                        .copied()
                        .map(|p| p as u32)
                        .filter(|&p| p < tail)
                        .collect();
                    sup.sort_unstable();
                    for_each_position_run(&sup, blocks, d, |b, first, len| {
                        block_work[b].accum.push(AccumRun {
                            term: ti as u32,
                            g0: first,
                            len,
                            // b_cls is a full-d buffer: class position `p`
                            // lives at slot `p`.
                            src: unsafe { b_data.add(first as usize) },
                            zero: false,
                        });
                    });
                    for (b, &(c0, c1)) in ranges.iter().enumerate() {
                        let lo = c0.max(tail);
                        if lo < c1 {
                            block_work[b].accum.push(AccumRun {
                                term: ti as u32,
                                g0: lo,
                                len: c1 - lo,
                                src: unsafe { b_data.add(lo as usize) },
                                zero: false,
                            });
                        }
                    }
                }
                TermSource::Node(k) => {
                    let ki = *k as usize - 2;
                    let (ptr, _) = compact_ptrs[ki];
                    // The node's exact batch scatter set (sorted class
                    // positions): the positions its sweep writes, hence
                    // the only slots that can be non-zero here.
                    let sup: Vec<u32> = self.nonzeros[ki].iter().map(|&p| p as u32).collect();
                    debug_assert!(sup.windows(2).all(|w| w[0] < w[1]));
                    let fused = node_terms[ki].len() == 1;
                    debug_assert!(!fused || node_terms[ki][0] == ti as u32);
                    for &CompactSlice { base, rs, re } in &layouts[ki] {
                        // The support positions inside the slice (both
                        // lists ascending); position `p` lives at slot
                        // `base + p - rs` and belongs to block
                        // `p * blocks / d`.
                        let s0 = sup.partition_point(|&p| p < rs);
                        let s1 = sup.partition_point(|&p| p < re);
                        for_each_position_run(&sup[s0..s1], blocks, d, |b, first, len| {
                            block_work[b].accum.push(AccumRun {
                                term: ti as u32,
                                g0: first,
                                len,
                                src: unsafe { ptr.add((base + first - rs) as usize) },
                                zero: fused,
                            });
                        });
                    }
                }
            }
        }
        // Zero runs for the slots no fused accum run zeroes: interior
        // nodes (never accumulated — shared brackets read only as sweep
        // operands) and, defensively, multi-rooted nodes.
        #[cfg(debug_assertions)]
        let mut zero_covered: Vec<u64> = vec![0; internal];
        for (ni, &(ptr, _)) in compact_ptrs.iter().enumerate() {
            if node_terms[ni].len() == 1 {
                continue;
            }
            let sup: Vec<u32> = self.nonzeros[ni].iter().map(|&p| p as u32).collect();
            debug_assert!(sup.windows(2).all(|w| w[0] < w[1]));
            for &CompactSlice { base, rs, re } in &layouts[ni] {
                let s0 = sup.partition_point(|&p| p < rs);
                let s1 = sup.partition_point(|&p| p < re);
                #[cfg(debug_assertions)]
                let covered = &mut zero_covered[ni];
                for_each_position_run(&sup[s0..s1], blocks, d, |b, first, len| {
                    #[cfg(debug_assertions)]
                    {
                        *covered += len as u64;
                    }
                    block_work[b].zero.push(ZeroRun {
                        ptr,
                        s0: base + first - rs,
                        len,
                    });
                });
            }
        }
        // Plan-time coverage invariant: every live slot of every node is
        // zeroed exactly once per fold — by its single referencing term's
        // fused runs (which cover the node's whole scatter set) or by the
        // explicit zero runs above.
        #[cfg(debug_assertions)]
        for (ni, sup) in self.nonzeros.iter().enumerate() {
            let fused = node_terms[ni].len() == 1;
            debug_assert_eq!(
                zero_covered[ni] as usize + if fused { sup.len() } else { 0 },
                sup.len(),
                "node {ni}: live-slot zero coverage incomplete"
            );
        }
        let work_shared = SendBlockWork(Arc::new(block_work));
        let mut block_closures: Vec<Box<dyn Fn(usize) + Send + Sync>> =
            Vec::with_capacity(2 * rhss.len() + 1);
        {
            let work = SendBlockWork(Arc::clone(&work_shared.0));
            let rhs_len = rhss[0].len();
            let rhs = SendConst(rhss[0].as_ptr());
            let b = SendRaw(b_data);
            let ranges = ranges.clone();
            block_closures.push(Box::new(move |bi: usize| {
                // Binding references first forces whole-value captures: a
                // direct `rhs.0.add(..)` place access would let the precise
                // capture split off the bare raw-pointer field, losing the
                // Send/Sync wrapper's promise.
                let (work, rhs, b, ranges) = (&work, &rhs, &b, &ranges);
                // zero this block's slot runs of the compact buffers
                for run in &work.0[bi].zero {
                    unsafe {
                        for si in run.s0 as usize..(run.s0 + run.len) as usize {
                            *run.ptr.add(si) = U::default();
                        }
                    }
                }
                let (c0, c1) = (ranges[bi].0 as usize, ranges[bi].1 as usize);
                #[allow(clippy::needless_range_loop)]
                unsafe {
                    for g in c0..c1 {
                        let pw = perm[g] as usize;
                        *b.0.add(g) = if pw < rhs_len {
                            (*rhs.0.add(pw)).clone()
                        } else {
                            U::default()
                        };
                    }
                }
            }));
        }
        for (fold_idx, rhs) in rhss.iter().enumerate() {
            {
                let rhs_len = rhs.len();
                let rhs = SendConst(rhs.as_ptr());
                let b = SendRaw(b_data);
                let ranges = ranges.clone();
                block_closures.push(Box::new(move |bi: usize| {
                    let (rhs, b) = (&rhs, &b);
                    let ranges = &ranges;
                    let (c0, c1) = (ranges[bi].0 as usize, ranges[bi].1 as usize);
                    #[allow(clippy::needless_range_loop)]
                    unsafe {
                        for g in c0..c1 {
                            let pw = perm[g] as usize;
                            *b.0.add(g) = if pw < rhs_len {
                                (*rhs.0.add(pw)).clone()
                            } else {
                                U::default()
                            };
                        }
                    }
                }));
            }
            {
                let structure = Arc::clone(&self.structure);
                let acc = SendRaw(acc_data);
                let b = SendRaw(b_data);
                let work = SendBlockWork(Arc::clone(&work_shared.0));
                let ranges = ranges.clone();
                // DEBUG: lets the last task snapshot acc/b_cls after all
                // ranges' accumulate+zero are complete.
                let dbg_ctr = Arc::new(std::sync::atomic::AtomicUsize::new(0));
                block_closures.push(Box::new(move |bi: usize| {
                    // Whole-value captures (see the lead closure above).
                    let (work, acc, b, structure, ranges) = (&work, &acc, &b, &structure, &ranges);
                    let terms = &structure.terms;
                    unsafe {
                        // This block's intersecting runs only, in the
                        // original term order: each position's add sequence
                        // matches the full walk bit for bit. Fused runs
                        // zero their source slot right after consuming it
                        // — the node's slots are dead once its single
                        // referencing term has read them — so the fold's
                        // zeroing touches the same cache lines the adds
                        // already do.
                        for run in &work.0[bi].accum {
                            let weight = &terms[run.term as usize].1;
                            let g0 = run.g0 as usize;
                            // SAFETY: the run's positions stay inside the
                            // block's disjoint range of `acc` and inside the
                            // source buffer's run of live slots; the fused
                            // zero-back writes only slots this run just read
                            // (owned by this block, dead afterwards).
                            if run.zero {
                                for (g, si) in
                                    (g0..g0 + run.len as usize).zip(0..run.len as usize)
                                {
                                    // `raw_mul`/`raw_add_assign_ptr` skip the
                                    // float wrappers' per-op NaN checks
                                    // (raw-float fast path, see lie-rs
                                    // `raw_mul`'s NaN policy).
                                    lie_rs::raw_add_assign_ptr(
                                        acc.0.add(g),
                                        &lie_rs::raw_mul(&*run.src.add(si), weight),
                                    );
                                    *run.src.add(si) = U::default();
                                }
                            } else {
                                for (g, si) in
                                    (g0..g0 + run.len as usize).zip(0..run.len as usize)
                                {
                                    lie_rs::raw_add_assign_ptr(
                                        acc.0.add(g),
                                        &lie_rs::raw_mul(&*run.src.add(si), weight),
                                    );
                                }
                            }
                        }
                        // zero this block's slot runs of the compact buffers
                        for run in &work.0[bi].zero {
                            for si in run.s0 as usize..(run.s0 + run.len) as usize {
                                *run.ptr.add(si) = U::default();
                            }
                        }
                    }
                    // DEBUG (BATCH_CKS): the last task of this stage
                    // snapshots acc+b_cls (all ranges' accumulate/zero/gather
                    // are ordered by the AcqRel counter RMWs; the snapshot
                    // runs before this task's publish, so the waiters — who
                    // proceed on the counter — cannot overtake it).
                    if lie_rs::CKS_ON.get().copied().unwrap_or(false)
                        && let Some(sink) = lie_rs::DEBUG_WRITES.get()
                    {
                        use std::sync::atomic::Ordering as O;
                        let n = dbg_ctr.fetch_add(1, O::AcqRel) + 1;
                        if n == ranges.len() {
                            let hh = |ptr: *const U, len: usize| {
                                let mut h: u64 = 0xcbf29ce484222325;
                                for k in 0..len {
                                    // SAFETY: debug read of the live
                                    // buffers, ordered by the AcqRel
                                    // counter above.
                                    let bits =
                                        unsafe { (*(ptr.add(k) as *const u64)).to_ne_bytes() };
                                    h ^= bits.iter().fold(0u64, |a, b| (a << 8) | *b as u64);
                                    h = h.wrapping_mul(0x100000001b3);
                                }
                                h
                            };
                            let fold = fold_idx;
                            let ha = hh(acc.0, ranges.last().map(|r| r.1 as usize).unwrap_or(0));
                            let hb = hh(b.0, ranges.last().map(|r| r.1 as usize).unwrap_or(0));
                            sink.lock()
                                .unwrap()
                                .push(format!("CKC fold={fold} acc={ha:016x} b={hb:016x}"));
                        }
                    }
                }));
            }
        }

        let mut block_stages: Vec<ClassBatchStage<U>> = Vec::with_capacity(2 * rhss.len() + 1);
        block_stages.push(ClassBatchStage::blocks(blocks, &*block_closures[0]));
        for f in 0..rhss.len() {
            block_stages.push(ClassBatchStage::blocks(blocks, &*block_closures[1 + 2 * f]));
            block_stages.push(ClassBatchStage::blocks(blocks, &*block_closures[2 + 2 * f]));
        }
        let mut stages: Vec<&ClassBatchStage<U>> =
            Vec::with_capacity(block_stages.len() + sweep_stages.len() * rhss.len());
        stages.push(&block_stages[0]);
        for f in 0..rhss.len() {
            stages.push(&block_stages[1 + 2 * f]);
            for s in &sweep_stages {
                stages.push(s);
            }
            stages.push(&block_stages[2 + 2 * f]);
        }

        // DEBUG (BATCH_CKS): publish the class-space working buffers for the
        // per-stage snapshots; cleared after the walk joins so a later
        // per-fold dispatch (which has no such buffers) never hashes stale
        // pointers.
        use std::sync::atomic::Ordering as DbgO;
        lie_rs::DEBUG_AB_ACC.store(acc_data as usize, DbgO::Relaxed);
        lie_rs::DEBUG_AB_B.store(b_data as usize, DbgO::Relaxed);
        lie_rs::DEBUG_AB_D.store(d, DbgO::Relaxed);
        // Fold units for the slot policy: the stage chain repeats its
        // gather/sweep/accumulate sub-chain once per displacement (plus
        // one leading zero stage), so each unit's barrier cost is paid
        // per fold — the policy normalizes by the unit count.
        run_class_batch(series, order, &stages, rhss.len().max(1));
        lie_rs::DEBUG_AB_ACC.store(0, DbgO::Relaxed);
        lie_rs::DEBUG_AB_B.store(0, DbgO::Relaxed);
        lie_rs::DEBUG_AB_D.store(0, DbgO::Relaxed);

        // Epilogue: the class-space accumulator back to public basis order
        // (assignment — it holds the full sum).
        let acc_cls = acc_cell.into_inner();
        for (k, &src) in inv.iter().enumerate().take(series.coefficients.len()) {
            series.coefficients[k] = acc_cls[src as usize].clone();
        }

        // The lists were used with (dense accumulator, common displacement
        // support) throughout; record those as the atom supports so the
        // next per-fold comparison is against reality.
        self.atom_a.clear();
        self.atom_a.extend_from_slice(a_nonzero);
        self.atom_b.clear();
        self.atom_b.extend_from_slice(b_nonzero);
        self.lists_built = true;
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
        for (source, weight) in &self.structure.terms {
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
        if self.dirty.len() != count || self.dirty.first().is_none_or(|w| w.len() != words) {
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
        *seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
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
}
