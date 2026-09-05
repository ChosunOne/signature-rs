//! SIMD-across-folds ("cohort") fast path for the commutation kernel.
//!
//! A cohort executes FOUR independent folds' sweeps as one walk of the
//! shared plan, carrying the folds' coefficient values as `[f64; 4]` (or
//! `[f32; 4]`) lane vectors over SoA-interleaved buffers: the four lanes'
//! value of buffer slot `s` lives at elements `[s*4 .. s*4+4]`, lane `l` at
//! `s*4 + l`. Decompositions, scatter indices, gating and BCH weights are
//! lane-invariant plan data (a fold plan is a pure function of the basis
//! and degree — the folds differ only in coefficient VALUES), so one walk
//! serves all four lanes and every position access is a single 4-lane
//! vector load/store.
//!
//! The lane arithmetic replicates the scalar kernel's per-lane operation
//! order exactly (plain mul + add, no FMA — Rust does not contract by
//! default), so cohort results are bit-identical to the scalar kernel's
//! per-fold results (0 ulp; verified by the spike and by the log-signature
//! fold's cohort-vs-scalar oracle test).
//!
//! The kernel is gated to the repr-transparent raw-float coefficient types
//! (`cohort_supported`): the lane vectors reinterpret the coefficient
//! slots through the primitive type exactly like the raw-float fast path
//! (`raw_ops`' NaN policy — overflowing inputs may produce NaN/Inf in
//! the slots; callers audit).
use super::*;
use ordered_float::NotNan;
use std::any::TypeId;

/// The four lanes of one cohort sweep.
pub const LANES: usize = 4;

/// True exactly for the coefficient types the cohort kernel supports:
/// `f32`/`f64` and their `NotNan` wrappers (repr-transparent over the
/// primitive — same set as the raw-float fast path). Every other
/// coefficient type keeps the scalar kernel.
pub fn supported<U: 'static>() -> bool {
    TypeId::of::<U>() == TypeId::of::<f64>()
        || TypeId::of::<U>() == TypeId::of::<NotNan<f64>>()
        || TypeId::of::<U>() == TypeId::of::<f32>()
        || TypeId::of::<U>() == TypeId::of::<NotNan<f32>>()
}

#[inline(always)]
fn is_f64_lane<U: 'static>() -> bool {
    TypeId::of::<U>() == TypeId::of::<f64>()
        || TypeId::of::<U>() == TypeId::of::<NotNan<f64>>()
}

#[inline(always)]
fn is_f32_lane<U: 'static>() -> bool {
    TypeId::of::<U>() == TypeId::of::<f32>()
        || TypeId::of::<U>() == TypeId::of::<NotNan<f32>>()
}

/// Reads the four lane values of one SoA buffer slot.
///
/// # Safety
/// `src` must be valid for [`LANES`] contiguous `U` elements, and `U`
/// must be a repr-transparent form of the lane primitive `L` — the
/// caller's TypeId dispatch guarantees both.
#[inline(always)]
unsafe fn load4<U, L>(src: *const U) -> [L; 4]
where
    L: Copy + Default + Add<Output = L> + Mul<Output = L> + Neg<Output = L> + 'static,
{
    // SAFETY: see the doc; the slots are 4-element-strided so an
    // unaligned read covers every layout the engine can produce.
    unsafe { (src as *const [L; LANES]).read_unaligned() }
}

/// Writes the four lane values of one SoA buffer slot.
///
/// # Safety
/// `dst` must be valid for [`LANES`] contiguous `U` elements, and `U`
/// must be a repr-transparent form of `L` (see [`load4`]). Written
/// values may be NaN for overflowing float inputs — the raw-float NaN
/// policy (callers audit).
#[inline(always)]
unsafe fn store4<U, L>(dst: *mut U, v: [L; LANES])
where
    L: Copy + Default + Add<Output = L> + Mul<Output = L> + Neg<Output = L> + 'static,
{
    // SAFETY: see the doc.
    unsafe { (dst as *mut [L; LANES]).write_unaligned(v) }
}

/// Broadcasts one coefficient value to all four lanes.
#[inline(always)]
fn splat4<U, L>(value: &U) -> [L; LANES]
where
    L: Copy + Default + Add<Output = L> + Mul<Output = L> + Neg<Output = L> + 'static,
{
    // SAFETY: a single live `U` element, repr-transparent over `L`.
    let s = unsafe { (value as *const U).cast::<L>().read() };
    [s; LANES]
}

#[inline(always)]
fn mul4<L: Copy + Mul<Output = L>>(a: [L; LANES], b: [L; LANES]) -> [L; LANES] {
    [
        a[0] * b[0],
        a[1] * b[1],
        a[2] * b[2],
        a[3] * b[3],
    ]
}

#[inline(always)]
fn add4<L: Copy + Add<Output = L>>(a: [L; LANES], b: [L; LANES]) -> [L; LANES] {
    [
        a[0] + b[0],
        a[1] + b[1],
        a[2] + b[2],
        a[3] + b[3],
    ]
}

#[inline(always)]
fn neg4<L: Copy + Neg<Output = L>>(a: [L; LANES]) -> [L; LANES] {
    [-a[0], -a[1], -a[2], -a[3]]
}

/// Cohort multiply-accumulate for the fold's accumulate phase:
/// `dst4 += weight × src4` for all four lanes, where `weight` is the
/// (shared) BCH weight of the term being accumulated and `src`/`dst`
/// are SoA slots. The per-lane operation order matches the scalar
/// engine's `raw_add_assign(dst, &raw_mul(src, weight))` exactly, so
/// cohort accumulation is bit-identical to the scalar kernel's.
///
/// # Safety
/// `dst`/`src` must each be valid for [`LANES`] contiguous `U`
/// elements, and the cohort capability ([`supported`]) must hold for
/// `U` (the caller's gate; other types panic).
pub unsafe fn add_mul4<U>(dst: *mut U, src: *const U, weight: &U)
where
    U: 'static,
{
    if is_f64_lane::<U>() {
        // SAFETY: the caller's contract; `U` is an f64 form.
        unsafe {
            let w = splat4::<U, f64>(weight);
            let v = load4::<U, f64>(src);
            let d = load4::<U, f64>(dst);
            store4::<U, f64>(dst, add4(d, mul4(v, w)));
        }
    } else if is_f32_lane::<U>() {
        // SAFETY: as the f64 branch, with f32 lanes.
        unsafe {
            let w = splat4::<U, f32>(weight);
            let v = load4::<U, f32>(src);
            let d = load4::<U, f32>(dst);
            store4::<U, f32>(dst, add4(d, mul4(v, w)));
        }
    } else {
        panic!(
            "cohort accumulate requires a raw repr-transparent float \
             coefficient type (f32/f64 or NotNan)"
        );
    }
}

/// The cohort (4-lane SoA) sweep of one pack: the exact
/// `FoldWalk::sweep_pack_range` single-phase per-unit visit order —
/// one term per active entry (computed inline, in table order) with
/// its row contributions added straight into the SoA result buffer —
/// with the four folds' values carried as lane vectors over
/// SoA-interleaved buffers (the jobs' operand/result pointers point at
/// SoA bases; lane `l` of slot `s` lives at element `s*LANES + l`).
///
/// Callers must have verified [`supported`] for `U`; the dispatch below
/// selects the lane width from the coefficient type and panics
/// otherwise (the cohort engine never forms for other types).
pub(super) fn sweep_pack_range<U>(
    entries: &[Entry],
    degrees: &[u8],
    coeffs: &[U],
    decomp_rels: &[u32],
    stage: &ClassBatchStage<'_, U>,
    pack: usize,
) where
    U: 'static,
{
    if is_f64_lane::<U>() {
        sweep_pack_impl::<U, f64>(entries, degrees, coeffs, decomp_rels, stage, pack);
    } else if is_f32_lane::<U>() {
        sweep_pack_impl::<U, f32>(entries, degrees, coeffs, decomp_rels, stage, pack);
    } else {
        panic!(
            "cohort sweep requires a raw repr-transparent float \
             coefficient type (f32/f64 or NotNan)"
        );
    }
}

/// The typed cohort sweep kernel. The visit order, unit/pack
/// structure, shift arithmetic and per-lane operation order replicate
/// `FoldWalk::sweep_pack_range` exactly; the only difference is that
/// every position access loads/stores [`LANES`] lanes at once.
fn sweep_pack_impl<U, L>(
    entries: &[Entry],
    degrees: &[u8],
    coeffs: &[U],
    decomp_rels: &[u32],
    stage: &ClassBatchStage<'_, U>,
    pack: usize,
) where
    L: Copy
        + Default
        + Add<Output = L>
        + Mul<Output = L>
        + Neg<Output = L>
        + Send
        + Sync
        + 'static,
    U: 'static,
{
    let (start, end) = stage.packs[pack];
    let jobs = stage.jobs;
    let mut t = start;
    while t < end {
        let ji = stage.tasks[t].0 as usize;
        let mut e = t + 1;
        while e < end && stage.tasks[e].0 == stage.tasks[t].0 {
            e += 1;
        }
        let job = &jobs[ji];
        let gating = &stage.gateways[ji];
        let a_shift = job.a_shift;
        let b_shift = job.b_shift;
        let r_shift = job.r_shift;
        for &(ji_task, ui) in &stage.tasks[t..e] {
            let _ = ji_task;
            let unit = &gating.units[ui as usize];
            let td = unit.td as usize;
            // SAFETY: the shift tables have >= max_degree+1 entries
            // and `td` is a degree of the table (see the plan). The
            // unit's degree-`td` slice starts at `rs` in the working
            // layout; the job's compact SoA buffer stores the slice at
            // `(rs - r_shift[td]) * LANES`.
            let base =
                (unit.rs as usize - unsafe { *r_shift.add(td) } as usize) * LANES;
            for tp in unit.ticket_start as usize..unit.ticket_end as usize {
                let ticket = gating.tickets[tp];
                let e_idx = (ticket & TICKET_INDEX_MASK) as usize;
                let entry = entries[e_idx];
                let p_active = ticket & TICKET_P_ACTIVE != 0;
                let q_active = ticket & TICKET_Q_ACTIVE != 0;
                let (i, j) = (entry.i as usize, entry.j as usize);
                let p = degrees[i] as usize;
                let q = degrees[j] as usize;
                // SAFETY: the shift tables have >= max_degree+1
                // entries and p/q are degrees of the table; i and j
                // are class positions in the operands' supports,
                // mapped into the SoA buffers at slots
                // `(position - shift) * LANES`.
                let term: [L; LANES] = unsafe {
                    if p_active {
                        let a_sh_p = *a_shift.add(p) as usize;
                        let b_sh_q = *b_shift.add(q) as usize;
                        let mut t4 = mul4(
                            load4::<U, L>(job.a.add((i - a_sh_p) * LANES)),
                            load4::<U, L>(job.b.add((j - b_sh_q) * LANES)),
                        );
                        if q_active {
                            let a_sh_q = *a_shift.add(q) as usize;
                            let b_sh_p = *b_shift.add(p) as usize;
                            t4 = add4(
                                t4,
                                neg4(mul4(
                                    load4::<U, L>(job.a.add((j - a_sh_q) * LANES)),
                                    load4::<U, L>(job.b.add((i - b_sh_p) * LANES)),
                                )),
                            );
                        }
                        t4
                    } else {
                        let a_sh_q = *a_shift.add(q) as usize;
                        let b_sh_p = *b_shift.add(p) as usize;
                        neg4(mul4(
                            load4::<U, L>(job.a.add((j - a_sh_q) * LANES)),
                            load4::<U, L>(job.b.add((i - b_sh_p) * LANES)),
                        ))
                    }
                };
                let from = entry.decomp_start as usize;
                let to = entries[e_idx + 1].decomp_start as usize;
                for (k, &rel) in decomp_rels[from..to].iter().enumerate() {
                    let c = &coeffs[from + k];
                    // SAFETY: single-writer 4-wide RMW (unit word-set
                    // ownership, packs never split a unit), in bounds
                    // (the gating addresses the job's compact SoA
                    // result space).
                    unsafe {
                        let d = load4::<U, L>(job.result.add(base + rel as usize * LANES));
                        store4::<U, L>(
                            job.result.add(base + rel as usize * LANES),
                            add4(d, mul4(splat4::<U, L>(c), term)),
                        );
                    }
                }
            }
        }
        t = e;
    }
}
