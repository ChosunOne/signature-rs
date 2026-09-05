//! The Arc copy-on-write node-list machinery: steady-list snapshots,
//! adoption from the plan cache and the level-parallel collecting rebuild.

use lie_rs::{ClassOrder, ClassOrderedCommutation, LieSeries};
use lyndon_rs::generators::Generator;
use num_traits::{One, Zero};
use std::hash::Hash;
use std::ops::{AddAssign, MulAssign, Neg};
use std::sync::Arc;

use super::{CommutatorDag, DagNode, node_slice};

/// Arc-shared immutable per-node scatter-target lists — the dag field
/// type and the leaf steady-plan cache's snapshot type. A dag adopts a
/// snapshot by `Arc` reference (zero copy); every mutation path is
/// copy-on-write (replace the `Arc` wholesale, or `Arc::unwrap_or_clone`
/// when exclusively held), so a shared snapshot is never mutated.
pub(crate) type SharedNodeLists = Arc<Vec<Vec<usize>>>;
/// Arc-shared immutable per-node dirty bitsets — shared exactly like
/// [`SharedNodeLists`].
pub(crate) type SharedDirtyLists = Arc<Vec<Vec<u64>>>;

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
    /// Snapshots the current steady-list state for plan caching: the node
    /// scatter lists and their dirty bitsets exactly as the last
    /// collecting rebuild left them (support-determined content — the
    /// pairing with the derived plan is the caller's, see
    /// [`Self::adopt_steady_lists`]). `Arc` clones: the snapshot shares
    /// this dag's lists until the dag first mutates them (copy-on-write).
    pub(crate) fn steady_lists_snapshot(
        &self,
    ) -> (SharedNodeLists, SharedDirtyLists) {
        (Arc::clone(&self.nonzeros), Arc::clone(&self.dirty))
    }

    /// Adopts a support-determined steady-list snapshot: the node scatter
    /// lists, dirty bitsets and atom supports produced by a previous
    /// [`Self::ensure_lists_steady`] fixed point for the SAME supports.
    /// Sound because collection is support-determined — the lists are a
    /// pure function of (table, class order, atom supports) — so the
    /// adopted state is bit-identical to what a fresh collecting rebuild
    /// on this dag would produce, minus the kernel sweeps. The snapshot's
    /// class order is installed too ([`Self::fold_batch`] requires a
    /// prior steady fold's order; [`Self::fold_batch_cohort`] would
    /// re-derive it).
    ///
    /// Validate-then-commit: returns `false` WITHOUT mutating `self` when
    /// the snapshot's shape does not fit this dag's structure (a pooled
    /// dag of a different configuration) or carries a different class
    /// order than the dag already holds — the caller falls back to a
    /// regular collecting rebuild. The snapshot's lists are installed by
    /// `Arc` reference — zero copy — and this dag's first list mutation
    /// thereafter copies (the shared snapshot is never mutated).
    pub(crate) fn adopt_steady_lists(
        &mut self,
        a_plan: &[usize],
        b_nonzero: &[usize],
        nonzeros: Arc<Vec<Vec<usize>>>,
        dirty: SharedDirtyLists,
        order: Arc<ClassOrder>,
    ) -> bool {
        if self.structure.nodes.len() - 2 != nonzeros.len() || dirty.len() != nonzeros.len() {
            return false;
        }
        if let Some(existing) = self.class_order.get()
            && !Arc::ptr_eq(existing, &order)
        {
            return false;
        }
        self.nonzeros = nonzeros;
        self.dirty = dirty;
        self.atom_a.clear();
        self.atom_a.extend_from_slice(a_plan);
        self.atom_b.clear();
        self.atom_b.extend_from_slice(b_nonzero);
        self.lists_built = true;
        let _ = self.class_order.set(order);
        true
    }

    /// Whether the node lists are the built fixed point for these atom
    /// supports — the list half of [`Self::batch_eligible`] without the
    /// accumulator-density requirement (a ONE-fold batch needs no density:
    /// its single sweep's gating is computed from the actual supports, and
    /// no later fold can outgrow it). This is the merge-cohort eligibility
    /// predicate; multi-fold cohort tails (leaf chunks) use the full
    /// [`Self::batch_eligible`].
    pub(crate) fn lists_steady_for(&self, a_nonzero: &[usize], b_nonzero: &[usize]) -> bool {
        self.lists_built && self.atom_a == a_nonzero && self.atom_b == b_nonzero
    }

    /// (Re)derives the node support lists when they are not the fixed point
    /// for these atom supports — the collecting pass of `evaluate`'s rebuild
    /// branch, run against THROWAWAY sweep results. Used by the merge cohort
    /// to bring a pooled DAG's plan up to the merge supports before the
    /// cohort sweep: the collected targets are value-independent (the
    /// kernel's scatter set is determined by the supports' gating alone), so
    /// one pass on one lane's values serves the whole cohort. The DAG's own
    /// node buffers are left with stale values — the cohort sweep uses its
    /// own compact buffers, and the next scalar `evaluate` re-zeroes them.
    pub(crate) fn ensure_lists_steady<T>(
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
        if self.lists_steady_for(a_nonzero, b_nonzero) {
            return;
        }
        let internal = self.structure.nodes.len() - 2;
        self.ensure_buffers(series.coefficients.len(), internal);
        let order = self.class_order.get_or_init(|| series.class_order());
        let inv = order.inv();
        let a_cls = series.class_coefficients(order, a);
        let b_cls = series.class_coefficients(order, b);
        let a_nz_cls: Vec<usize> = a_nonzero.iter().copied().map(|i| inv[i] as usize).collect();
        let b_nz_cls: Vec<usize> = b_nonzero.iter().copied().map(|i| inv[i] as usize).collect();
        // Per-level PARALLEL collecting sweep. The rebuild is the fold
        // pipeline's serial mountain otherwise: one node's collect kernel
        // call with dense supports costs O(active entries), the DAG has as
        // many nodes as BCH terms, and every plan derivation and
        // support-change rebuild walks all of them node-by-node (measured
        // ~28ms serial at 4x8 — larger than the parallel batch dispatch it
        // feeds). A level's nodes depend only on strictly earlier levels'
        // lists and write disjoint (scratch, dirty, list) triples, so each
        // level is one conflict-free parallel batch — the same level
        // structure the batch dispatch uses. Collection is
        // support-determined: the schedule cannot change the lists.
        let d = series.coefficients.len();
        // Copy-on-write: a dag sharing an adopted snapshot's lists must
        // never mutate them through the `Arc` — the level-parallel tasks
        // below write through raw pointers into these locals. Unwrap when
        // exclusively held (free), clone when shared.
        let mut lists = Arc::unwrap_or_clone(std::mem::take(&mut self.nonzeros));
        let mut dirties = Arc::unwrap_or_clone(std::mem::take(&mut self.dirty));
        {
            let nodes = &self.structure.nodes;
            let buffers = &self.buffers;
            // The collect tasks capture ONLY Send+Sync data: the T-free
            // kernel inputs (feasible table + class order + max degree +
            // basis length — all basis-derived integer tables), the U
            // operand buffers, and the per-node slot pointers. The letter
            // type never crosses a thread (same contract as the batch
            // dispatch's T-free walk). Each node appears in exactly one
            // level and is collected exactly once per call, so the per-node
            // slots are written by exactly one task; a task reads only its
            // nodes' children's slots (strictly earlier levels — already
            // final). This mirrors the SeriesTemplate safety contract in
            // signature-rs, minus the T exposure: no T-bearing value is
            // reachable from the tasks at all.
            struct NodeSlots {
                lists: Vec<*mut Vec<usize>>,
                dirty: Vec<*mut Vec<u64>>,
            }
            unsafe impl Send for NodeSlots {}
            unsafe impl Sync for NodeSlots {}
            let slots = NodeSlots {
                lists: lists.iter_mut().map(|l| l as *mut Vec<usize>).collect(),
                dirty: dirties.iter_mut().map(|x| x as *mut Vec<u64>).collect(),
            };
            // Bind through a shared reference so the task closures capture
            // `&NodeSlots` (Sync via the unsafe impl above) instead of the
            // per-field `&Vec<*mut _>` borrows, whose element raw pointers
            // are not Sync.
            let slots = &slots;
            let table = series.feasible_decompositions_handle().clone();
            let basis_len = series.basis.len();
            let max_degree = series.max_degree;
            use rayon::prelude::*;
            for level in self.structure.levels.iter().skip(1) {
                level.par_iter().for_each(|&k| {
                    let (left, right) = match nodes[k as usize] {
                        DagNode::Binary { left, right } => (left, right),
                        DagNode::Atom(_) => unreachable!("atoms are the first two nodes"),
                    };
                    let lbuf = node_slice(left, buffers, &a_cls, &b_cls);
                    let rbuf = node_slice(right, buffers, &a_cls, &b_cls);
                    let child_list = |id: u32| -> &[usize] {
                        if id < 2 {
                            if id == 0 { &a_nz_cls } else { &b_nz_cls }
                        } else {
                            // SAFETY: the child's slot was published by an
                            // earlier level (the level ordering is the
                            // serial walk's topological order) and is not
                            // written again this call.
                            unsafe { &*slots.lists[id as usize - 2] }
                        }
                    };
                    let lnz = child_list(left);
                    let rnz = child_list(right);
                    // SAFETY: this node's slot is written only here (one
                    // task per node per call — levels partition the nodes).
                    // The raw pointers are copied to locals first: `&mut *`
                    // through the shared `slots` place would demand a
                    // mutable borrow of it.
                    let list_ptr = slots.lists[k as usize - 2];
                    let list = unsafe { &mut *list_ptr };
                    let dirty_ptr = slots.dirty[k as usize - 2];
                    let dirty = unsafe { &mut *dirty_ptr };
                    list.clear();
                    // One scratch result buffer per node call: the sweep's
                    // value outputs are discarded; only the collected
                    // scatter-target list is kept.
                    let mut scratch = vec![U::default(); d];
                    LieSeries::<T, U>::class_collect_kernel(
                        &table,
                        basis_len,
                        max_degree,
                        order,
                        lbuf,
                        lnz,
                        rbuf,
                        rnz,
                        &mut scratch,
                        dirty,
                        list,
                    );
                });
            }
        }
        self.nonzeros = Arc::new(lists);
        self.dirty = Arc::new(dirties);
        self.atom_a.clear();
        self.atom_a.extend_from_slice(a_nonzero);
        self.atom_b.clear();
        self.atom_b.extend_from_slice(b_nonzero);
        self.lists_built = true;
    }
}
