//! Compiled evaluation plan for the BCH concatenation expression: shared
//! subtrees interned to node ids, evaluated as a straight-line loop over
//! reused per-node buffers instead of recursive per-fold evaluation.
//!
//! Module map:
//! - [`structure`]: the node graph (`DagNode`, `DagStructure`, `TermSource`)
//!   and the `from_bch_series` construction, including `strip_dead_nodes`.
//! - [`block_work`]: per-block work lists and raw-pointer Send wrappers
//!   shared by both batch engines.
//! - [`evaluate`]: single-fold evaluation, batch eligibility and the
//!   term-accumulation epilogue.
//! - [`fold_batch`]: the scalar batch engine.
//! - [`fold_cohort`]: the 4-lane SIMD-across-folds cohort engine, its lane
//!   type, kill switch and counters.
//! - [`steady_lists`]: the Arc copy-on-write node-list machinery.
//! - [`pool`]: `clone_shallow` plan sharing for the dag pool.
//! - [`tests`]: the module's test suite.

mod block_work;
mod evaluate;
mod fold_batch;
mod fold_cohort;
mod pool;
mod steady_lists;
mod structure;
#[cfg(test)]
mod tests;

pub(crate) use fold_cohort::{CohortLane, cohort_capable, cohort_type_capable};
#[cfg(test)]
pub(crate) use fold_cohort::{COHORT_ENGINE_RUNS, COHORT_LANE_FOLDS, set_cohort_off};
pub(crate) use steady_lists::{SharedDirtyLists, SharedNodeLists};
pub(crate) use structure::{DagNode, DagStructure, TermSource};

use lie_rs::{ClassOrder, GatingCache};
use num_traits::{One, Zero};
use std::hash::Hash;
use std::ops::{AddAssign, MulAssign, Neg};
use std::sync::{Arc, OnceLock};

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
    ///
    /// `Arc`-shared: adoption from the leaf steady-plan cache installs the
    /// cache's immutable snapshot by reference (zero copy), and every
    /// mutation path is copy-on-write — the owning `Arc` is replaced
    /// wholesale (rebuilds/re-records) or unwrapped when exclusively held
    /// (`Arc::unwrap_or_clone`) — so a shared snapshot is never mutated
    /// through the `Arc`.
    pub(crate) nonzeros: SharedNodeLists,
    /// Deduplication scratch for the collecting rebuild, one bitset per
    /// internal node sized to the basis. Shared/copied exactly like
    /// [`Self::nonzeros`].
    dirty: SharedDirtyLists,
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
            self.nonzeros = Arc::new((0..count).map(|_| Vec::new()).collect::<Vec<_>>());
        }
        let words = d.div_ceil(64);
        if self.dirty.len() != count || self.dirty.first().is_none_or(|w| w.len() != words) {
            self.dirty = Arc::new((0..count).map(|_| vec![0u64; words]).collect::<Vec<_>>());
        }
    }
}
