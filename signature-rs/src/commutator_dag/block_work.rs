//! Per-block work lists and raw-pointer Send wrappers shared by the
//! scalar and cohort batch engines' block stages.

use std::sync::Arc;

/// Raw-pointer wrapper for the batch's cross-stage working buffers. All
/// access goes through these pointers (never through a `&mut` taken after
/// derivation); disjointness comes from the block ranges within a stage,
/// and the batch's stage counters order every cross-stage dependency.
pub(super) struct SendRaw<T>(pub(super) *mut T);
unsafe impl<T> Send for SendRaw<T> {}
unsafe impl<T> Sync for SendRaw<T> {}

/// The per-fold displacement `(pointer, length)` table for the batch's
/// gather task: the task reads fold `f`'s slice through it (raw-pointer
/// bearing, same Send/Sync argument as [`SendRaw`]: the pointers
/// target the caller's displacement storage, read-only for the whole
/// dispatch).
pub(super) struct SendRhsTable<U>(pub(super) Vec<(*const U, usize)>);
unsafe impl<U> Send for SendRhsTable<U> {}
unsafe impl<U> Sync for SendRhsTable<U> {}

/// One (cohort step, lane) displacement view for the cohort gather task:
/// lane `l`'s `f`-th remaining displacement slice, or null when the lane
/// has folded past step `f` (masked). Read-only for the whole dispatch —
/// same Send/Sync argument as [`SendRhsTable`].
pub(super) struct RhsView<U> {
    pub(super) ptr: *const U,
    pub(super) len: usize,
}

/// The cohort gather task's per-(step, lane) displacement views (raw-
/// pointer bearing, same Send/Sync argument as [`SendRhsTable`]).
pub(super) struct SendRhsViews<U>(pub(super) Arc<Vec<[RhsView<U>; 4]>>);
unsafe impl<U> Send for SendRhsViews<U> {}
unsafe impl<U> Sync for SendRhsViews<U> {}

/// One accumulate run of a block's work list: a contiguous (class
/// position, source slot) stretch of one BCH term whose positions fall in
/// the block's range.
pub(super) struct AccumRun<U> {
    /// Index into `DagStructure::terms` for the run's BCH weight.
    pub(super) term: u32,
    /// First class position of the run.
    pub(super) g0: u32,
    /// Number of class positions in the run.
    pub(super) len: u32,
    /// Source pointer: the node's compact slot holding position `g0` (the
    /// fold displacement at `g0` for the displacement term); advances with
    /// `g0`. Mutable because fused runs write `U::default()` back after
    /// consuming each slot (the pointer is derived from the batch's
    /// mutable compact/displacement buffers).
    pub(super) src: *mut U,
    /// Fused zero-after-add: each source slot is written back to
    /// `U::default()` as it is consumed. Set only for a node referenced by
    /// exactly one term (its slots are dead once that term has read them);
    /// the displacement term's `b_cls` source is rewritten wholesale by
    /// the next fold's gather, so it never needs this.
    pub(super) zero: bool,
}

/// One zero run of a block's work list: a contiguous stretch of one
/// compact node buffer's slots inside the block's class range.
pub(super) struct ZeroRun<U> {
    /// Target compact buffer pointer.
    pub(super) ptr: *mut U,
    /// First compact slot of the run.
    pub(super) s0: u32,
    /// Number of compact slots in the run.
    pub(super) len: u32,
}

/// A block's per-batch work lists: the accumulate runs and zero runs
/// intersecting the block's class-position range (see `fold_batch`).
pub(super) struct BlockWork<U> {
    pub(super) accum: Vec<AccumRun<U>>,
    pub(super) zero: Vec<ZeroRun<U>>,
}

/// The per-block work lists of one batch, shared by the batch's block
/// closures through an `Arc` (raw-pointer bearing, same Send/Sync argument
/// as [`SendBufPtrs`]: the pointers target stable frame-owned allocations
/// whose cross-stage ordering the stage counters provide).
pub(super) struct SendBlockWork<U>(pub(super) std::sync::Arc<Vec<BlockWork<U>>>);
unsafe impl<U> Send for SendBlockWork<U> {}
unsafe impl<U> Sync for SendBlockWork<U> {}

/// Groups a sorted ascending position list into maximal runs of consecutive
/// positions within one block (position `p` belongs to block
/// `p * blocks / d`), invoking `emit(block, first_position, run_length)`
/// per run in ascending order. Consecutive positions of one block merge
/// into a single run; a block boundary or a position gap starts a new one.
pub(super) fn for_each_position_run(
    positions: &[u32],
    step: usize,
    mut emit: impl FnMut(usize, u32, u32),
) {
    let mut i = 0;
    while i < positions.len() {
        let first = positions[i];
        let block = first as usize / step;
        let mut last = first;
        let mut j = i + 1;
        while j < positions.len()
            && positions[j] == last + 1
            && positions[j] as usize / step == block
        {
            last = positions[j];
            j += 1;
        }
        emit(block, first, last - first + 1);
        i = j;
    }
}
