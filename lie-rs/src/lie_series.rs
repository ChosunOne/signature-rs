use std::any::TypeId;
use std::collections::{HashMap, HashSet};
use std::fmt::{Debug, Display};
use std::hash::Hash;
use std::marker::PhantomData;
use std::ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign};
use std::sync::atomic::{AtomicU32, AtomicUsize, Ordering};
use std::sync::{Arc, OnceLock};

use commutator_rs::{Commutator, CommutatorTerm, comm};
use lyndon_rs::generators::Generator;
use lyndon_rs::lyndon::LyndonWord;
use num_traits::{One, Zero};

use crate::feasible_decompositions::{
    self, Entry, FeasibleDecompositions, UnitRange, TICKET_INDEX_MASK, TICKET_P_ACTIVE,
    TICKET_Q_ACTIVE,
};

// Re-exported: the class-contiguous ordering handle behind
// `ClassOrderedCommutation`.
pub use crate::feasible_decompositions::ClassOrder;

pub use self::raw_ops::{raw_add_assign, raw_add_assign_ptr, raw_mul};

/// Raw-float fast path for the commutation kernel's coefficient arithmetic.
///
/// For the float coefficient types used in production the kernel's hot loops
/// run their `*`/`+=` through these helpers instead of the arithmetic impls
/// of the checking [`ordered_float::NotNan`] wrapper (each of which pays a
/// per-operation NaN check, `ucomisd`+`jp` in the sweep's inner loop). The
/// wrapper is `#[repr(transparent)]` over the primitive, so the helpers
/// reinterpret the operands and compute in the primitive type; the
/// `TypeId` comparisons constant-fold per monomorphization, so the dispatch
/// is free for concrete callers and every other coefficient type (e.g.
/// `Ratio<i128>`) keeps the identical generic path.
///
/// # NaN policy
/// For the raw-float instantiations no per-operation NaN check runs. Finite
/// inputs cannot produce NaN through `*`/`+=`; overflow produces infinities
/// exactly as the checked path does, and only `inf - inf`-style combinations
/// of those can yield NaN. Such a NaN is stored bit-for-bit into the
/// `NotNan` slot (an invalid value for the wrapper's logical invariant, but
/// not undefined behavior — the wrapper has no niche). Callers that persist
/// results in `NotNan` slots audit them once per fold step (the
/// log-signature fold) and fail loudly instead of leaving the broken
/// invariant behind. Results are otherwise bitwise identical to the checked
/// path.
mod raw_ops {
    use super::*;
    use ordered_float::NotNan;

    /// `a * b` without the float wrappers' per-operation NaN checks (see the
    /// module-level [NaN policy]).
    #[inline(always)]
    pub fn raw_mul<U>(a: &U, b: &U) -> U
    where
        U: Clone + Mul<Output = U> + 'static,
    {
        if TypeId::of::<U>() == TypeId::of::<NotNan<f64>>() {
            // SAFETY: `U` is `NotNan<f64>` (check above), which is
            // `#[repr(transparent)]` over `f64`; every `NotNan<f64>` is a
            // valid `f64`, so reading the operands through `f64` pointers is
            // sound. The result may be NaN for overflowing inputs — the
            // documented NaN policy covers that (callers audit).
            unsafe {
                let r = *(a as *const U).cast::<f64>() * *(b as *const U).cast::<f64>();
                std::ptr::read((&r as *const f64).cast::<U>())
            }
        } else if TypeId::of::<U>() == TypeId::of::<NotNan<f32>>() {
            // SAFETY: as the `f64` branch, with `NotNan<f32>` over `f32`.
            unsafe {
                let r = *(a as *const U).cast::<f32>() * *(b as *const U).cast::<f32>();
                std::ptr::read((&r as *const f32).cast::<U>())
            }
        } else {
            a.clone() * b.clone()
        }
    }

    /// `*dst += src` without the float wrappers' per-operation NaN checks
    /// (see the module-level [NaN policy]).
    #[inline(always)]
    pub fn raw_add_assign<U>(dst: &mut U, src: &U)
    where
        U: Clone + AddAssign + 'static,
    {
        // SAFETY: `dst` is a live, uniquely borrowed `U`.
        unsafe { raw_add_assign_ptr(dst as *mut U, src as *const U) }
    }

    /// Raw-pointer counterpart of [`raw_add_assign`]: the parallel sweeps'
    /// scatter targets are accessed through raw pointers (disjoint write
    /// regions across tasks).
    ///
    /// # Safety
    /// `dst` must be valid for writes and point to a live `U`; `src` must be
    /// valid for reads and point to a live `U`. Values written may be NaN for
    /// overflowing float inputs (the documented NaN policy — callers audit).
    #[inline(always)]
    pub unsafe fn raw_add_assign_ptr<U>(dst: *mut U, src: *const U)
    where
        U: Clone + AddAssign + 'static,
    {
        if TypeId::of::<U>() == TypeId::of::<NotNan<f64>>() {
            // SAFETY: `U` is `NotNan<f64>` (check above), `#[repr(transparent)]`
            // over `f64`; both pointers reference live values and `dst` is
            // uniquely owned by the caller's contract. A NaN write (overflow
            // cancellation) is covered by the NaN policy.
            unsafe {
                *(dst.cast::<f64>()) += *src.cast::<f64>();
            }
        } else if TypeId::of::<U>() == TypeId::of::<NotNan<f32>>() {
            // SAFETY: as the `f64` branch, with `NotNan<f32>` over `f32`.
            unsafe {
                *(dst.cast::<f32>()) += *src.cast::<f32>();
            }
        } else {
            // SAFETY: `dst`/`src` are live per the caller's contract; the
            // generic path is the wrapper's own checked `+=`.
            unsafe { (*dst) += (*src).clone() };
        }
    }
}

/// SIMD-across-folds ("cohort") fast path for the commutation kernel.
///
/// A cohort executes FOUR independent folds' sweeps as one walk of the
/// shared plan, carrying the folds' coefficient values as `[f64; 4]` (or
/// `[f32; 4]`) lane vectors over SoA-interleaved buffers: the four lanes'
/// value of buffer slot `s` lives at elements `[s*4 .. s*4+4]`, lane `l` at
/// `s*4 + l`. Decompositions, scatter indices, gating and BCH weights are
/// lane-invariant plan data (a fold plan is a pure function of the basis
/// and degree — the folds differ only in coefficient VALUES), so one walk
/// serves all four lanes and every position access is a single 4-lane
/// vector load/store.
///
/// The lane arithmetic replicates the scalar kernel's per-lane operation
/// order exactly (plain mul + add, no FMA — Rust does not contract by
/// default), so cohort results are bit-identical to the scalar kernel's
/// per-fold results (0 ulp; verified by the spike and by the log-signature
/// fold's cohort-vs-scalar oracle test).
///
/// The kernel is gated to the repr-transparent raw-float coefficient types
/// (`cohort_supported`): the lane vectors reinterpret the coefficient
/// slots through the primitive type exactly like the raw-float fast path
/// (`raw_ops`' NaN policy — overflowing inputs may produce NaN/Inf in
/// the slots; callers audit).
pub mod cohort {
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
}

#[cfg(feature = "progress")]
use indicatif::{ProgressBar, ProgressStyle};

/// Represents a formal power series in the free Lie algebra.
///
/// A `LieSeries` is a linear combination of Lyndon words (represented as commutator terms)
/// with coefficients from a ring. This structure provides the foundation for computations
/// involving Baker-Campbell-Hausdorff series and other Lie algebraic operations.
pub struct LieSeries<T, U> {
    /// The Lyndon word basis for the free Lie algebra.
    pub basis: Arc<Vec<LyndonWord<T>>>,
    /// The commutator representation of the Lyndon basis elements.
    pub commutator_basis: Arc<Vec<CommutatorTerm<U, T>>>,
    /// Maps arbitrary commutator terms to their decomposition in the basis.
    /// Invariant: no stored decomposition has a zero coefficient.
    pub(crate) feasible_decompositions: Arc<FeasibleDecompositions<U>>,

    /// Maps basis elements to their index positions for efficient lookup.
    pub commutator_basis_index_map: Arc<HashMap<u64, usize>>,
    /// The coefficients corresponding to each basis element.
    pub coefficients: Vec<U>,
    /// The maximum degree of terms included in this series.
    pub max_degree: usize,
}

impl<T: Debug, U: Debug> Debug for LieSeries<T, U> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("LieSeries")
            .field("basis", &self.basis)
            .field("commutator_basis", &self.commutator_basis)
            .field(
                "feasible_decompositions",
                &self.feasible_decompositions.len(),
            )
            .field(
                "commutator_basis_index_map",
                &self.commutator_basis_index_map,
            )
            .field("coefficients", &self.coefficients)
            .field("max_degree", &self.max_degree)
            .finish()
    }
}

impl<T: Display, U: Display + One + PartialEq> Display for LieSeries<T, U> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for (i, (coefficient, basis_term)) in self
            .coefficients
            .iter()
            .zip(self.commutator_basis.iter())
            .enumerate()
        {
            if i == 0 {
                write!(f, "{coefficient} {basis_term}")?;
                continue;
            }
            write!(f, " + {coefficient} {basis_term}")?;
        }
        Ok(())
    }
}

impl<T: Clone, U: Clone> Clone for LieSeries<T, U> {
    fn clone(&self) -> Self {
        Self {
            basis: self.basis.clone(),
            commutator_basis: Arc::clone(&self.commutator_basis),
            commutator_basis_index_map: Arc::clone(&self.commutator_basis_index_map),
            coefficients: self.coefficients.clone(),
            feasible_decompositions: Arc::clone(&self.feasible_decompositions),
            max_degree: self.max_degree,
        }
    }
}

impl<T, U> Index<usize> for LieSeries<T, U> {
    type Output = U;

    fn index(&self, index: usize) -> &Self::Output {
        &self.coefficients[index]
    }
}

impl<
    T: Clone + Ord + Generator + Hash,
    U: Clone + One + Zero + Eq + MulAssign + Neg<Output = U> + Hash,
> Index<LyndonWord<T>> for LieSeries<T, U>
{
    type Output = U;

    fn index(&self, index: LyndonWord<T>) -> &Self::Output {
        let term = CommutatorTerm::<U, T>::from(&index);
        let i = self.commutator_basis_index_map[&term.unit_hash()];
        &self.coefficients[i]
    }
}

impl<
    T: Clone + Ord + Generator + Hash,
    U: Clone + One + Zero + Eq + MulAssign + Neg<Output = U> + Hash,
> Index<&LyndonWord<T>> for LieSeries<T, U>
{
    type Output = U;

    fn index(&self, index: &LyndonWord<T>) -> &Self::Output {
        let term = CommutatorTerm::<U, T>::from(index);
        let i = self.commutator_basis_index_map[&term.unit_hash()];
        &self.coefficients[i]
    }
}

impl<T, U> IndexMut<usize> for LieSeries<T, U> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.coefficients[index]
    }
}

impl<
    T: Clone + Ord + Generator + Hash,
    U: Clone + One + Zero + Eq + MulAssign + Neg<Output = U> + Hash,
> IndexMut<LyndonWord<T>> for LieSeries<T, U>
{
    fn index_mut(&mut self, index: LyndonWord<T>) -> &mut Self::Output {
        let term = CommutatorTerm::<U, T>::from(&index);
        let i = self.commutator_basis_index_map[&term.unit_hash()];
        &mut self.coefficients[i]
    }
}

impl<
    T: Clone + Ord + Generator + Hash,
    U: Clone + One + Zero + Eq + MulAssign + Neg<Output = U> + Hash,
> IndexMut<&LyndonWord<T>> for LieSeries<T, U>
{
    fn index_mut(&mut self, index: &LyndonWord<T>) -> &mut Self::Output {
        let term = CommutatorTerm::<U, T>::from(index);
        let i = self.commutator_basis_index_map[&term.unit_hash()];
        &mut self.coefficients[i]
    }
}

impl<T: Clone, U: Clone + Add<Output = U>> Add for &LieSeries<T, U> {
    type Output = LieSeries<T, U>;

    fn add(self, rhs: Self) -> Self::Output {
        let mut result = self.clone();
        for i in 0..self.coefficients.len() {
            result.coefficients[i] = self.coefficients[i].clone() + rhs.coefficients[i].clone();
        }
        result
    }
}

impl<T, U: Clone + Add<Output = U>> Add for LieSeries<T, U> {
    type Output = LieSeries<T, U>;

    fn add(mut self, rhs: Self) -> Self::Output {
        for i in 0..self.coefficients.len() {
            self.coefficients[i] = self.coefficients[i].clone() + rhs.coefficients[i].clone();
        }
        self
    }
}

impl<T, U: Clone + AddAssign> AddAssign for LieSeries<T, U> {
    fn add_assign(&mut self, rhs: Self) {
        for i in 0..self.coefficients.len() {
            self.coefficients[i] += rhs.coefficients[i].clone();
        }
    }
}

impl<T, U: Clone + AddAssign> AddAssign<&Self> for LieSeries<T, U> {
    fn add_assign(&mut self, rhs: &Self) {
        for i in 0..self.coefficients.len() {
            self.coefficients[i] += rhs.coefficients[i].clone();
        }
    }
}

impl<T: Clone, U: Clone + Sub<Output = U>> Sub for &LieSeries<T, U> {
    type Output = LieSeries<T, U>;

    fn sub(self, rhs: Self) -> Self::Output {
        let mut result = self.clone();
        for i in 0..self.coefficients.len() {
            result.coefficients[i] = self.coefficients[i].clone() - rhs.coefficients[i].clone();
        }
        result
    }
}

impl<T, U: Clone + Sub<Output = U>> Sub for LieSeries<T, U> {
    type Output = LieSeries<T, U>;

    fn sub(mut self, rhs: Self) -> Self::Output {
        for i in 0..self.coefficients.len() {
            self.coefficients[i] = self.coefficients[i].clone() - rhs.coefficients[i].clone();
        }
        self
    }
}

impl<T, U: Clone + SubAssign> SubAssign for LieSeries<T, U> {
    fn sub_assign(&mut self, rhs: Self) {
        for i in 0..self.coefficients.len() {
            self.coefficients[i] -= rhs.coefficients[i].clone();
        }
    }
}

impl<T, U: Clone + SubAssign> SubAssign<&Self> for LieSeries<T, U> {
    fn sub_assign(&mut self, rhs: &Self) {
        for i in 0..self.coefficients.len() {
            self.coefficients[i] -= rhs.coefficients[i].clone();
        }
    }
}

impl<T, U: Clone + Mul<Output = U> + MulAssign> Mul<U> for LieSeries<T, U> {
    type Output = Self;

    fn mul(mut self, rhs: U) -> Self::Output {
        for i in 0..self.coefficients.len() {
            self.coefficients[i] *= rhs.clone();
        }
        self
    }
}

impl<T: Clone, U: Clone + Mul<Output = U> + MulAssign> Mul<U> for &LieSeries<T, U> {
    type Output = LieSeries<T, U>;

    fn mul(self, rhs: U) -> Self::Output {
        let mut result = self.clone();
        for i in 0..self.coefficients.len() {
            result.coefficients[i] *= rhs.clone();
        }
        result
    }
}

impl<T: Clone, U: Clone + Mul<Output = U> + MulAssign> MulAssign<U> for LieSeries<T, U> {
    fn mul_assign(&mut self, rhs: U) {
        for c in &mut self.coefficients {
            *c *= rhs.clone();
        }
    }
}

/// Cheap deterministic per-basis setup shared by every construction path:
/// the commutator-basis terms, their unit-hash membership set, and the
/// unit-hash -> basis-index map.
struct CommutatorPrologue<U, T> {
    commutator_basis: Vec<CommutatorTerm<U, T>>,
    commutator_basis_set: HashSet<u64>,
    commutator_basis_index_map: HashMap<u64, usize>,
}

impl<
    T: Clone + Ord + Generator + Hash + Eq,
    U: Clone + One + Zero + Eq + MulAssign + Neg<Output = U> + Hash + AddAssign + Ord,
> LieSeries<T, U>
{
    #[must_use]
    #[cfg_attr(
        feature = "tracing",
        tracing::instrument(skip_all, fields(basis_len = basis.len()))
    )]
    /// DEBUG/telemetry: per-position Lyndon degree and the degree-slice
    /// starts of the feasible-decomposition table (used by support-size
    /// probes in downstream crates).
    #[doc(hidden)]
    pub fn debug_degree_layout(&self) -> (Vec<u8>, Vec<u32>) {
        self.feasible_decompositions.debug_degree_layout()
    }

    /// Cheap deterministic per-basis setup shared by every construction
    /// path: the commutator-basis terms, their unit-hash membership set,
    /// and the unit-hash -> basis-index map.
    fn commutator_prologue(basis: &[LyndonWord<T>]) -> CommutatorPrologue<U, T> {
        let mut commutator_basis = Vec::<CommutatorTerm<U, T>>::with_capacity(basis.len());
        // Only the tracing build reports it; kept for the debug! shorthand.
        #[cfg_attr(not(feature = "tracing"), allow(unused_variables))]
        let max_degree = if basis.is_empty() {
            0
        } else {
            basis[basis.len() - 1].len()
        };
        for word in basis {
            commutator_basis.push(CommutatorTerm::from(word));
        }
        #[cfg(feature = "tracing")]
        tracing::debug!(
            max_degree,
            basis_len = commutator_basis.len(),
            "converted Lyndon basis to commutator terms"
        );
        let commutator_basis_set = commutator_basis
            .iter()
            .map(commutator_rs::CommutatorTerm::unit_hash)
            .collect::<HashSet<_>>();

        let mut commutator_basis_index_map = HashMap::new();
        for (i, a) in commutator_basis.iter().enumerate() {
            commutator_basis_index_map.insert(a.unit_hash(), i);
        }
        #[cfg(feature = "tracing")]
        tracing::debug!("built commutator basis maps");
        CommutatorPrologue {
            commutator_basis,
            commutator_basis_set,
            commutator_basis_index_map,
        }
    }

    /// Builds the feasible-decomposition (structure-constant) table for a
    /// basis: for every feasible pair `(i, j)` the commutator `[w_i, w_j]`
    /// is rewritten into the Lyndon basis. Deterministic — pairs are
    /// visited in index order and each decomposition's output terms are
    /// sorted — so the table is a pure function of the basis (plus the
    /// coefficient type's arithmetic). This is the expensive part of
    /// [`LieSeries::new`]; it is factored out so callers can build the
    /// table once and share it via
    /// [`LieSeries::with_feasible_decompositions`].
    #[must_use]
    pub fn build_feasible_decompositions(basis: &[LyndonWord<T>]) -> FeasibleDecompositions<U> {
        let prologue = Self::commutator_prologue(basis);
        Self::structure_constants(
            basis,
            &prologue.commutator_basis,
            &prologue.commutator_basis_set,
            &prologue.commutator_basis_index_map,
        )
    }

    /// The structure-constant pass over the pair table. `basis` must be the
    /// basis `commutator_basis` / `commutator_basis_set` /
    /// `commutator_basis_index_map` were derived from.
    fn structure_constants(
        basis: &[LyndonWord<T>],
        commutator_basis: &[CommutatorTerm<U, T>],
        commutator_basis_set: &HashSet<u64>,
        commutator_basis_index_map: &HashMap<u64, usize>,
    ) -> FeasibleDecompositions<U> {
        let max_degree = if basis.is_empty() {
            0
        } else {
            basis[basis.len() - 1].len()
        };
        // Per-word letter contents: the table is grouped by (target degree,
        // content class), which requires each word's letter multiset.
        let mut alphabet: Vec<T> = basis
            .iter()
            .filter(|w| w.len() == 1)
            .map(|w| w.letters[0].clone())
            .collect();
        alphabet.sort();
        alphabet.dedup();
        let basis_contents: Vec<Vec<u8>> = basis
            .iter()
            .map(|w| {
                let mut c = vec![0u8; alphabet.len()];
                for l in &w.letters {
                    let k = alphabet
                        .binary_search(l)
                        .expect("basis letter missing from alphabet");
                    c[k] += 1;
                }
                c
            })
            .collect();
        let mut feasible_builder = feasible_decompositions::Builder::new(&basis_contents);

        #[cfg(feature = "tracing")]
        let _structure_constants_span = tracing::debug_span!(
            "compute_structure_constants",
            pairs = basis.len() * basis.len(),
            max_degree
        )
        .entered();
        #[cfg(feature = "progress")]
        let pb = {
            let style = ProgressStyle::with_template(
                "[{elapsed_precise}] [{bar:35.cyan/blue}] {pos}/{len} structure-constant rows ({msg})",
            )
            .unwrap()
            .progress_chars("=>-");
            let pb = ProgressBar::new(basis.len() as u64).with_style(style);
            pb
        };
        for (i, a) in commutator_basis.iter().enumerate() {
            #[cfg(feature = "progress")]
            {
                pb.set_message(format!("degree {}", a.degree()));
                pb.inc(1);
            }
            #[cfg(feature = "tracing")]
            if i == 0 || commutator_basis[i - 1].degree() != a.degree() {
                tracing::debug!(
                    row = i,
                    degree = a.degree(),
                    "computing structure constants"
                );
            }
            for (j, b) in commutator_basis.iter().enumerate() {
                if j <= i || max_degree < a.degree() + b.degree() {
                    continue;
                }
                #[cfg(feature = "tracing")]
                tracing::trace!(i, j, "computing commutator term");
                let mut term = comm![a, b];
                #[cfg(feature = "tracing")]
                tracing::trace!(i, j, "sorting commutator term");
                term.lyndon_sort();

                // For some non-basis term T, and its decomposition to basis terms A, B, and C, ...
                // T -> [c_1 * A, c_2 * B, c_3 * C, ...]
                #[cfg(feature = "tracing")]
                tracing::trace!(i, j, "decomposing into Lyndon basis");
                let basis_terms = term
                    .lyndon_basis_decomposition(&commutator_basis_set)
                    .into_iter()
                    .filter(|x| !x.coefficient().is_zero())
                    .collect::<Vec<_>>();
                // Get the coefficients [c_1, c_2, c_3, ...]
                let basis_term_coefficients = basis_terms
                    .iter()
                    .map(|x| x.coefficient().clone())
                    .collect::<Vec<_>>();
                // Get the indices i of A, B, C, ...
                // [i_A, i_B, i_C, ...]
                let basis_term_indices = basis_terms
                    .into_iter()
                    .map(|x| commutator_basis_index_map[&x.unit_hash()])
                    .collect::<Vec<_>>();
                #[cfg(feature = "tracing")]
                tracing::trace!(
                    i,
                    j,
                    terms = basis_term_indices.len(),
                    "decomposed commutator into basis terms"
                );
                feasible_builder.push(i, j, &basis_term_indices, &basis_term_coefficients);
            }
        }

        #[cfg(feature = "tracing")]
        tracing::debug!(
            basis_len = basis.len(),
            "finished computing structure constants"
        );
        #[cfg(feature = "progress")]
        pb.finish_with_message("done");
        feasible_builder.finish()
    }

    /// Builds an empty series over `basis` with a pre-built
    /// feasible-decomposition table, skipping the structure-constant
    /// computation. The table must have been produced for exactly this
    /// basis (see [`LieSeries::build_feasible_decompositions`]); sharing
    /// one table across series over the same basis is sound because the
    /// table is read-only after construction.
    pub fn with_feasible_decompositions(
        basis: Arc<Vec<LyndonWord<T>>>,
        coefficients: Vec<U>,
        feasible_decompositions: Arc<FeasibleDecompositions<U>>,
    ) -> Self {
        let prologue = Self::commutator_prologue(&basis);
        let max_degree = if basis.is_empty() {
            0
        } else {
            basis[basis.len() - 1].len()
        };
        Self {
            basis,
            commutator_basis: Arc::new(prologue.commutator_basis),
            commutator_basis_index_map: Arc::new(prologue.commutator_basis_index_map),
            coefficients,
            feasible_decompositions,
            max_degree,
        }
    }

    /// Identity of the shared structure-constant table (diagnostics and
    /// tests: series over the same cached plan report the same id).
    #[doc(hidden)]
    #[must_use]
    pub fn table_id(&self) -> usize {
        Arc::as_ptr(&self.feasible_decompositions) as usize
    }

    pub fn new(basis: Vec<LyndonWord<T>>, coefficients: Vec<U>) -> Self {
        let prologue = Self::commutator_prologue(&basis);
        let feasible_decompositions = Self::structure_constants(
            &basis,
            &prologue.commutator_basis,
            &prologue.commutator_basis_set,
            &prologue.commutator_basis_index_map,
        );
        let max_degree = if basis.is_empty() {
            0
        } else {
            basis[basis.len() - 1].len()
        };
        Self {
            basis: Arc::new(basis),
            commutator_basis: Arc::new(prologue.commutator_basis),
            commutator_basis_index_map: Arc::new(prologue.commutator_basis_index_map),
            coefficients,
            feasible_decompositions: Arc::new(feasible_decompositions),
            max_degree,
        }
    }
}

impl<T, U> LieSeries<T, U> {
    #[must_use]
    pub fn with_coefficients(&self, coefficients: Vec<U>) -> Self {
        Self {
            basis: Arc::clone(&self.basis),
            commutator_basis: Arc::clone(&self.commutator_basis),
            feasible_decompositions: Arc::clone(&self.feasible_decompositions),
            commutator_basis_index_map: Arc::clone(&self.commutator_basis_index_map),
            coefficients,
            max_degree: self.max_degree,
        }
    }

    /// Decomposition of the bracket `[basis[i], basis[j]]` onto the basis,
    /// for pairs stored in the feasible table (`i != j` and degree-feasible).
    /// Returns parallel slices (basis indices, non-zero coefficients).
    #[must_use]
    pub fn decomposition(&self, i: usize, j: usize) -> Option<(&[u32], &[U], bool)> {
        self.feasible_decompositions.get(i, j)
    }
}

/// Target size of a parallel bundle, in visited entries.
const BUNDLE_TARGET_ENTRIES: usize = 2048;
/// Floor for the thread-adaptive bundle target, in visited entries: below
/// this a bundle's work is smaller than its dispatch cost, and unit
/// integrity (bundles never split a unit) provides the natural minimum.
const MIN_BUNDLE_ENTRIES: usize = 16;

/// The kernel prologue's output: the active units of the feasible-
/// decomposition table — each holding its active-entry ticket range, its
/// exact word set (the result positions it writes), and its contribution
/// (pair) count — plus the total pair count (the planner's work unit).
///
/// The division of parallel work is per unit: the unit word sets
/// PARTITION the sweep's write slots (units are unique by (target degree,
/// content) and every entry producing word `w` lives in `w`'s one unit),
/// so any pack cut between units is race-free, and a word's contributions
/// accumulate inside its one unit in table-entry order — exactly the
/// serial sweep's per-word float summation order, independent of pack
/// cuts, slot counts and claim interleaving.
///
/// Cloning is cheap (all hot fields are `Arc` slices) — the cohort sweep
/// stages clone their planned stage's gateways.
#[derive(Clone)]
struct KernelGating {
    /// Active units in table order, each with its ticket range, word-set
    /// range (into [`Self::unit_words`]) and pair weight. Units whose
    /// active entries all have empty rows are omitted (no work).
    units: Arc<[UnitRange]>,
    /// The flat active-entry ticket list, in table order (a subsequence of
    /// the entry stream — the entry table's natural locality). Each unit's
    /// `ticket_start .. ticket_end` slices this list; a unit's sweep
    /// streams its slice sequentially, computing each term once and adding
    /// its row contributions straight into the result buffer.
    tickets: Arc<[u32]>,
    /// The flat ascending active-word position list: the union of the
    /// units' word sets, sorted (units of one degree are ordered by
    /// content bytes, not position, so the union must be sorted for the
    /// sink/scatter-set walk order to stay byte-identical to the per-word
    /// fan-in's). Exactly the set of positions the sweep writes — the
    /// batch scatter sets are this list.
    unit_words: Arc<[u32]>,
    /// Total active contribution (pair) count: the sweep's work unit (one
    /// multiply-add each), used for pack balancing and the slot policy.
    total_pairs: usize,
}

/// Memo for the commutation kernel's prologues: maps a job's operand
/// support lists to its resolved gating — the active per-word write classes
/// with their `(ticket, decomp position)` contribution lists, presence
/// results baked in.
///
/// The key is a 128-bit fingerprint per side of the *exact* support lists
/// (not just their degree masks): the cached tickets bake per-entry
/// presence results, which degree masks do not determine — two supports
/// can cover the same degrees through different words (exactly what the
/// old per-entry bit tests re-resolved on every kernel call). The
/// fingerprints mix the list lengths in and avalanche-finish both streams,
/// so distinct lists reuse an entry only with negligible probability; an
/// entry's gating is reused solely for identical supports.
///
/// In a log-signature fold the DAG's node support lists are a value-
/// independent fixed point, so steady-state folds (and whole batches)
/// hit the cache for every job: the per-entry resolution runs once per
/// distinct (job, operand supports), not per fold.
///
/// The cache is valid only for the decomposition table of the series whose
/// prologues populated it.
#[derive(Clone, Default)]
pub struct GatingCache {
    /// Open-addressed linear-probe table keyed by the two support
    /// fingerprints. Distinct pairs per configuration are few (the DAG's
    /// distinct node-support pairs), so a fixed-capacity table with a
    /// cheap multiplicative hash beats a full hash map: the lookup runs per
    /// kernel call, and at small grids the SipHash + bucket walk was a
    /// measurable share of the fold.
    slots: Vec<Slot>,
}

/// Key + value for one [`GatingCache`] slot. `None` key = empty; the all-
/// zero fingerprint pair is a legitimate key (two empty supports), so
/// emptiness is tracked out of band.
#[derive(Clone, Default)]
struct Slot {
    key: Option<([u64; 2], [u64; 2])>,
    value: (Arc<[UnitRange]>, Arc<[u32]>, Arc<[u32]>, usize),
}

impl GatingCache {
    /// Looks up the support fingerprint pair, returning the memoized
    /// `(active units, flat tickets, active word positions, pair count)`.
    #[inline]
    fn get(
        &self,
        key: ([u64; 2], [u64; 2]),
    ) -> Option<&(Arc<[UnitRange]>, Arc<[u32]>, Arc<[u32]>, usize)> {
        let cap = self.slots.len();
        if cap == 0 {
            return None;
        }
        let mut idx = Self::hash(&key) & (cap - 1);
        loop {
            let slot = &self.slots[idx];
            match &slot.key {
                Some(k) if *k == key => return Some(&slot.value),
                Some(_) => idx = (idx + 1) & (cap - 1),
                None => return None,
            }
        }
    }

    /// Inserts, growing the table to the next power of two when full.
    fn insert(
        &mut self,
        key: ([u64; 2], [u64; 2]),
        value: (Arc<[UnitRange]>, Arc<[u32]>, Arc<[u32]>, usize),
    ) {
        if self.slots.len() * 3 / 4 <= self.len() {
            let old = std::mem::take(&mut self.slots);
            self.slots = vec![Slot::default(); (old.len() * 2).max(8)];
            for slot in old {
                if let Some(k) = slot.key {
                    self.insert_unchecked(k, slot.value);
                }
            }
        }
        self.insert_unchecked(key, value);
    }

    fn insert_unchecked(
        &mut self,
        key: ([u64; 2], [u64; 2]),
        value: (Arc<[UnitRange]>, Arc<[u32]>, Arc<[u32]>, usize),
    ) {
        let cap = self.slots.len();
        let mut idx = Self::hash(&key) & (cap - 1);
        while self.slots[idx].key.is_some() {
            idx = (idx + 1) & (cap - 1);
        }
        self.slots[idx] = Slot {
            key: Some(key),
            value,
        };
    }

    fn len(&self) -> usize {
        self.slots.iter().filter(|s| s.key.is_some()).count()
    }

    /// Multiplicative fold of the four mask words — cheaper than SipHash
    /// and collision-free enough for a table holding a handful of entries.
    #[inline]
    fn hash(key: &([u64; 2], [u64; 2])) -> usize {
        let mut h = 0xcbf29ce484222325u64;
        for w in [&key.0, &key.1] {
            for x in w {
                h = (h ^ x).wrapping_mul(0x100000001b3);
            }
        }
        h as usize
    }
}

/// Shared handle for the parallel sweep's result writes.
///
/// SAFETY (of the `Send`/`Sync` impls): bundles of units own disjoint
/// word sets — a basis word's (degree, content) determines the single
/// unit whose entries write it, and bundles never split a unit — so
/// concurrent read-modify-writes through `ptr` touch disjoint addresses.
/// `U: Send` covers the values moving across threads.
struct RawResult<'a, U> {
    ptr: *mut U,
    _marker: PhantomData<&'a mut [U]>,
}

unsafe impl<U: Send> Send for RawResult<'_, U> {}
unsafe impl<U: Send> Sync for RawResult<'_, U> {}

/// Chooses the entry table for a sweep. Monomorphized per layout, so each
/// path's inner loop compiles to exactly one code shape.
trait ScatterLayout {
    fn entries<'a>(class: Option<&'a ClassOrder>, direct: &'a [Entry]) -> &'a [Entry];
}

/// Public basis order: public entries, writes straight into the job's
/// result buffer.
struct DirectLayout;

/// Fully class-contiguous order: class-ordered operands (read through the
/// relabeled entry table) and a class-ordered result buffer — no scratch,
/// no permutation.
struct ClassInternalLayout;

impl ScatterLayout for DirectLayout {
    #[inline(always)]
    fn entries<'a>(_class: Option<&'a ClassOrder>, direct: &'a [Entry]) -> &'a [Entry] {
        direct
    }
}

impl ScatterLayout for ClassInternalLayout {
    #[inline(always)]
    fn entries<'a>(class: Option<&'a ClassOrder>, _direct: &'a [Entry]) -> &'a [Entry] {
        class
            .expect("class layout without a class order")
            .entries_cls()
    }
}

/// Parallel single-phase per-unit sweep: each bundle's units are visited
/// independently — one term per active entry, its row contributions added
/// straight into the result buffer (no intermediate term buffers, no phase
/// barrier). Units are atomic and their word sets partition the write
/// slots, so concurrent bundles never touch the same word; each word's
/// adds happen inside its one unit in table-entry order, which is the
/// per-word float summation sequence the serial sweep produces, so the
/// bits are unchanged. `L` selects the entry table layout — public
/// order into the job's result buffer, or class-contiguous order into a
/// class-ordered result buffer.
fn sweep_bundles_parallel<L: ScatterLayout, T, U>(
    jobs: &[KernelJob<'_, U>],
    writers: &[RawResult<U>],
    gateways: &[KernelGating],
    entries_slice: &[Entry],
    decomp_rels: &[u32],
    class_order: Option<&ClassOrder>,
    decomp_coeffs_slice: &[U],
    unit_bundles: &[Vec<(usize, u32)>],
) where
    T: Clone + Ord + Generator + Hash,
    U: Clone
        + Default
        + One
        + Zero
        + Eq
        + MulAssign
        + Neg<Output = U>
        + Hash
        + Mul<Output = U>
        + AddAssign
        + Send
        + Sync
        + 'static,
{
    use rayon::prelude::*;

    #[cfg(feature = "tracing")]
    let sweep_span = tracing::debug_span!(
        "kernel_sweep_parallel",
        jobs = jobs.len(),
        bundles = unit_bundles.len(),
        threads = rayon::current_num_threads(),
    )
    .entered();
    let entries_tbl = <L as ScatterLayout>::entries(class_order, entries_slice);

    // Single-phase per-unit sweep: each bundle's units are visited
    // independently — one term per active entry (computed inline, exactly
    // the old phase-1 formula) and its row contributions added straight
    // into the result buffer. Units are atomic and their word sets
    // partition the write slots — bundles never split one, so concurrent
    // bundles never touch the same word; each word's adds happen inside
    // its one unit in table-entry order, which is exactly the old
    // two-phase sweep's per-word float summation sequence (and the serial
    // sweep's), so the bits are unchanged. Accumulate-into semantics are
    // preserved (the buffer's current value starts every sum).
    unit_bundles
        .par_iter()
        .enumerate()
        .for_each(|(_bundle_index, bundle)| {
            #[cfg(feature = "tracing")]
            let _bundle_span = tracing::debug_span!(
                "kernel_bundle",
                bundle = _bundle_index,
                units = bundle.len(),
            )
            .entered();
            for &(ji, ui) in bundle {
                let job = &jobs[ji];
                let writer = &writers[ji];
                let gating = &gateways[ji];
                let unit = &gating.units[ui as usize];
                let base = unit.rs as usize;
                for tp in unit.ticket_start as usize..unit.ticket_end as usize {
                    let ticket = gating.tickets[tp];
                    let e = (ticket & TICKET_INDEX_MASK) as usize;
                    let entry = entries_tbl[e];
                    let p_active = ticket & TICKET_P_ACTIVE != 0;
                    let q_active = ticket & TICKET_Q_ACTIVE != 0;
                    let (i, j) = (entry.i as usize, entry.j as usize);
                    // SAFETY: i and j are positions < the operand lengths
                    // (the tickets resolved presence against the operand
                    // supports). `raw_mul` skips the float wrappers'
                    // per-op NaN checks (raw-float fast path); `-` never
                    // checks.
                    let term = unsafe {
                        if p_active {
                            let mut t = raw_mul(&*job.a.add(i), &*job.b.add(j));
                            if q_active {
                                raw_add_assign(
                                    &mut t,
                                    &-raw_mul(&*job.a.add(j), &*job.b.add(i)),
                                );
                            }
                            t
                        } else {
                            -raw_mul(&*job.a.add(j), &*job.b.add(i))
                        }
                    };
                    let from = entry.decomp_start as usize;
                    let to = entries_tbl[e + 1].decomp_start as usize;
                    for (k, &rel) in decomp_rels[from..to].iter().enumerate() {
                        // SAFETY: the unit's word set partitions the
                        // result positions (the ownership invariant) and
                        // this bundle owns the unit whole — single-writer
                        // RMW, in bounds (the gating addresses the working
                        // layout's result space).
                        unsafe {
                            raw_add_assign(
                                &mut *writer.ptr.add(base + rel as usize),
                                &raw_mul(&decomp_coeffs_slice[from + k], &term),
                            );
                        }
                    }
                }
            }
        });
    #[cfg(feature = "tracing")]
    drop(sweep_span);
}

impl<U> RawResult<'_, U>
where
    U: Clone + AddAssign + 'static,
{
}

/// One independent commutation: operand slices plus the destination buffer.
///
/// SAFETY contract for `LieSeries::commutator_coefficients_batch`: the
/// `result` buffers of the jobs passed to one batch call must be pairwise
/// disjoint (in a fold these are distinct DAG-node buffers). Within a job,
/// the anagram partition makes the units conflict-free.
pub struct KernelJob<'a, U> {
    /// The left operand's coefficients, as a raw pointer because the
    /// DAG-fold batch mutates the operand buffers through `UnsafeCell`
    /// between stages while these slices stay live — a shared reference
    /// would be frozen over concurrently-mutated memory (undefined
    /// behavior, and it lets the compiler cache operand reads across
    /// stage boundaries). Valid for `a_len` elements for the duration of
    /// the call.
    pub a: *const U,
    /// The left operand's length.
    pub a_len: usize,
    // SAFETY: `Sync` is sound because the batch's tasks only *read* the
    // operand buffers (through the raw pointers, ordered by the stage
    // counters) and write through `result` at indices owned by exactly
    // one task (disjoint buffers across jobs, anagram-disjoint units
    // within a job); `U: Send` covers the written values crossing threads.
    //
    // (The unsafe impls follow the struct definition.)
    /// The left operand's non-zero indices (superset of its support).
    pub a_nonzero: &'a [usize],
    /// The right operand's coefficients (see `a`).
    pub b: *const U,
    /// The right operand's length.
    pub b_len: usize,
    /// The right operand's non-zero indices.
    pub b_nonzero: &'a [usize],
    /// The destination buffer, as a raw pointer because the batch's parallel
    /// tasks write disjoint regions of it concurrently. Must remain valid
    /// for `result_len` elements for the duration of the call.
    pub result: *mut U,
    /// The destination buffer's length.
    pub result_len: usize,
    /// Compact-layout address shifts, indexed by Lyndon degree (tables must
    /// have at least `max_degree + 1` entries). Class position `x` of degree
    /// `d` lives at `x - shift[d]` in the corresponding buffer; a full-d
    /// buffer uses the all-zero `IDENTITY_SHIFTS` table. The batch
    /// fold's per-node buffers store only the degree slices the node's
    /// sweep writes (4-6x smaller than full-d for deep DAGs), so every
    /// operand/result access subtracts the hoisted per-degree shift.
    pub a_shift: *const u32,
    pub b_shift: *const u32,
    pub r_shift: *const u32,
}

/// The identity shift table for full-d buffers (see [`KernelJob`]).
pub const IDENTITY_SHIFTS: [u32; 128] = [0u32; 128];

// SAFETY: see the comment inside `KernelJob` — tasks write only indices
// owned by their job, and only values of `U: Send`.
unsafe impl<U: Send> Send for KernelJob<'_, U> {}
unsafe impl<U: Send> Sync for KernelJob<'_, U> {}

/// Batch kernel: evaluates several independent commutation jobs in one
/// parallel dispatch. Jobs must have disjoint `result` buffers; per-word
/// accumulation order matches the serial sweep exactly.
pub fn commutator_coefficients_batch<T, U>(
    a_series: &LieSeries<T, U>,
    jobs: &mut [KernelJob<'_, U>],
) where
    T: Clone + Ord + Generator + Hash,
    U: Clone
        + Default
        + One
        + Zero
        + Eq
        + MulAssign
        + Neg<Output = U>
        + Hash
        + Mul<Output = U>
        + AddAssign
        + Send
        + Sync
        + 'static,
{
    let mut cache = GatingCache::default();
    commutator_coefficients_batch_with_cache(a_series, jobs, &mut cache);
}

/// [`commutator_coefficients_batch`] with a caller-owned [`GatingCache`]
/// that persists across calls, amortizing the gating scan over repeated
/// jobs with equal degree support.
pub fn commutator_coefficients_batch_with_cache<T, U>(
    a_series: &LieSeries<T, U>,
    jobs: &mut [KernelJob<'_, U>],
    cache: &mut GatingCache,
) where
    T: Clone + Ord + Generator + Hash,
    U: Clone
        + Default
        + One
        + Zero
        + Eq
        + MulAssign
        + Neg<Output = U>
        + Hash
        + Mul<Output = U>
        + AddAssign
        + Send
        + Sync
        + 'static,
{
    // Prologue per job (serial): presence bitsets, degree masks, active
    // units (memoized). Cheap relative to the sweep.
    #[cfg(feature = "tracing")]
    let prologue_span = tracing::debug_span!("kernel_prologue", jobs = jobs.len()).entered();
    let gateways: Vec<KernelGating> = jobs
        .iter()
        .map(|j| {
            LieSeries::<T, U>::kernel_prologue_cached(a_series, j.a_nonzero, j.b_nonzero, cache)
        })
        .collect();
    let total: usize = gateways.iter().map(|g| g.total_pairs).sum();
    #[cfg(feature = "tracing")]
    {
        drop(prologue_span);
        tracing::debug!(
            jobs = jobs.len(),
            total_pairs = total,
            "kernel_prologue done"
        );
    }

    // Hoisted flat-table views: every contribution's ticket indexes the
    // entry table directly, so neither sweep touches the unit table again.
    let table = &a_series.feasible_decompositions;
    let entries_slice = table.entries();
    let decomp_coeffs_slice = table.decomp_coeffs();

    // With more than one thread available the batch always dispatches to
    // the parallel sweep; rayon's work stealing balances the pieces. Only a
    // single-threaded pool runs the serial sweep.
    let threads = rayon::current_num_threads();
    if threads == 1 {
        #[cfg(feature = "tracing")]
        let sweep_span = tracing::debug_span!(
            "kernel_sweep_serial",
            jobs = jobs.len(),
            total_pairs = total,
            threads
        )
        .entered();
        for (job, gating) in jobs.iter_mut().zip(&gateways) {
            // SAFETY: the job's result buffer is valid for `result_len`
            // elements (the struct's contract) and is exclusively ours here.
            let result = unsafe { std::slice::from_raw_parts_mut(job.result, job.result_len) };
            // SAFETY: these serial paths are exclusive per job (no
            // concurrent mutation), so shared slices over the raw operand
            // pointers are sound here.
            let (a_slice, b_slice) = unsafe {
                (
                    std::slice::from_raw_parts(job.a, job.a_len),
                    std::slice::from_raw_parts(job.b, job.b_len),
                )
            };
            // The gating's word classes carry public positions: one
            // single-writer accumulation + store per word, straight into
            // the public result buffer. (The old class-contiguous scratch
            // densified a scatter of repeated `+=`s; one store per word
            // needs no densifying, and the scratch + permutation epilogue
            // it paid is gone.)
            LieSeries::sweep_words_serial(
                a_series,
                a_slice,
                b_slice,
                gating,
                result,
                entries_slice,
                table.decomp_indices_rel(),
                &mut NoSink,
            );
        }
        #[cfg(feature = "tracing")]
        drop(sweep_span);
        return;
    }

    // SAFETY: each job's sweep target is a distinct buffer (disjoint across
    // jobs by the caller's contract), and within a job each word class is
    // written by exactly one task (bundles never split a class), so the
    // concurrent single-writer stores never alias. The buffers are not
    // otherwise accessed during the sweep.
    let writers: Vec<RawResult<U>> = jobs
        .iter()
        .map(|j| RawResult {
            ptr: j.result,
            _marker: PhantomData,
        })
        .collect();

    // Flatten (job, unit) pairs into bundles of roughly
    // `BUNDLE_TARGET_ENTRIES` contributions, weighted by the unit's pair
    // count. Units stay whole within a bundle (target/degree order,
    // preserving each word's accumulation context); a unit is never split.
    let mut bundles: Vec<Vec<(usize, u32)>> = vec![Vec::new()];
    let mut cur = 0usize;
    // Enough bundles that every thread can hold a piece, without dropping
    // below the per-task break-even size. The cut happens between units —
    // always safe: unit word sets never overlap.
    //
    // Balancing by decomposition volume (the alternative: bundles of equal
    // summed decomposition length) was measured and rejected in the
    // unit-atomic engine: per-bundle sweep time is ~uniform per work unit
    // across degree mixes, so counts are already the right balance metric.
    let bundle_target = (total / (2 * threads)).clamp(MIN_BUNDLE_ENTRIES, BUNDLE_TARGET_ENTRIES);
    for (ji, gating) in gateways.iter().enumerate() {
        for (ui, unit) in gating.units.iter().enumerate() {
            let pairs = unit.pairs as usize;
            if cur >= bundle_target {
                bundles.push(Vec::new());
                cur = 0;
            }
            bundles.last_mut().unwrap().push((ji, ui as u32));
            cur += pairs;
        }
    }

    sweep_bundles_parallel::<DirectLayout, T, U>(
        jobs,
        &writers,
        &gateways,
        entries_slice,
        table.decomp_indices_rel(),
        None,
        decomp_coeffs_slice,
        &bundles,
    )
}

/// Batch kernel in the class-contiguous working mode: every job's operand
/// slices are class-ordered, its support lists class-indexed, and its
/// result buffer receives the class-ordered sum directly — no scratch, no
/// permutation. Requires `order` to describe `a_series`' basis.
pub fn commutator_coefficients_class_batch_with_cache<T, U>(
    a_series: &LieSeries<T, U>,
    order: &ClassOrder,
    jobs: &mut [KernelJob<'_, U>],
    cache: &mut GatingCache,
) where
    T: Clone + Ord + Generator + Hash,
    U: Clone
        + Default
        + One
        + Zero
        + Eq
        + MulAssign
        + Neg<Output = U>
        + Hash
        + Mul<Output = U>
        + AddAssign
        + Send
        + Sync
        + 'static,
{
    debug_assert_eq!(
        order.inv().len(),
        a_series.basis.len(),
        "class ordering does not describe this series' basis"
    );

    #[cfg(feature = "tracing")]
    let prologue_span = tracing::debug_span!("kernel_prologue_class", jobs = jobs.len()).entered();
    let gateways: Vec<KernelGating> = jobs
        .iter()
        .map(|j| {
            LieSeries::<T, U>::kernel_prologue_cached_class(
                a_series,
                j.a_nonzero,
                j.b_nonzero,
                order,
                cache,
            )
        })
        .collect();
    let total: usize = gateways.iter().map(|g| g.total_pairs).sum();
    #[cfg(feature = "tracing")]
    drop(prologue_span);

    let table = &a_series.feasible_decompositions;
    let entries_slice = table.entries();
    let decomp_coeffs_slice = table.decomp_coeffs();
    let co: &ClassOrder = order;

    // Single-threaded pools sweep serially; the result buffer is already
    // class-ordered, so the sweep writes it directly.
    let threads = rayon::current_num_threads();
    if threads == 1 {
        for (job, gating) in jobs.iter_mut().zip(&gateways) {
            // SAFETY: the job's result buffer is valid for `result_len`
            // elements (the struct's contract) and is exclusively ours here;
            // the serial path is exclusive per job, so shared slices over
            // the raw operand pointers are sound.
            let result = unsafe { std::slice::from_raw_parts_mut(job.result, job.result_len) };
            let (a_slice, b_slice) = unsafe {
                (
                    std::slice::from_raw_parts(job.a, job.a_len),
                    std::slice::from_raw_parts(job.b, job.b_len),
                )
            };
            LieSeries::sweep_words_serial(
                a_series,
                a_slice,
                b_slice,
                gating,
                result,
                co.entries_cls(),
                co.decomp_cls(),
                &mut NoSink,
            );
        }
        return;
    }

    // SAFETY: each job's result is a distinct buffer (disjoint across jobs
    // by the caller's contract), and within a job each word class is
    // written by exactly one task (bundles never split a class), so the
    // concurrent single-writer stores never alias.
    let writers: Vec<RawResult<U>> = jobs
        .iter()
        .map(|j| RawResult {
            ptr: j.result,
            _marker: PhantomData,
        })
        .collect();

    // Flatten (job, unit) pairs into bundles — same cut rule as the
    // public kernel: pair-count balanced, cuts only between units (always
    // safe: unit word sets never overlap).
    let mut bundles: Vec<Vec<(usize, u32)>> = vec![Vec::new()];
    let mut cur = 0usize;
    let bundle_target = (total / (2 * threads)).clamp(MIN_BUNDLE_ENTRIES, BUNDLE_TARGET_ENTRIES);
    for (ji, gating) in gateways.iter().enumerate() {
        for (ui, unit) in gating.units.iter().enumerate() {
            let pairs = unit.pairs as usize;
            if cur >= bundle_target {
                bundles.push(Vec::new());
                cur = 0;
            }
            bundles.last_mut().unwrap().push((ji, ui as u32));
            cur += pairs;
        }
    }

    // The class-ordered result buffer needs no epilogue: the sweep's writes
    // are already final.
    sweep_bundles_parallel::<ClassInternalLayout, T, U>(
        jobs,
        &writers,
        &gateways,
        entries_slice,
        co.decomp_cls(),
        Some(co),
        decomp_coeffs_slice,
        &bundles,
    );
}

/// Class-contiguous working mode for the commutation kernel.
///
/// The class ordering depends only on the basis: every series over the
/// same basis rearranges coefficients identically, so one [`ClassOrder`]
/// handle — cheap to clone, `Arc` internally — amortizes the O(basis)
/// planning across operand series, across every kernel call of a fold,
/// and across batches of folds. Obtain the handle once via
/// [`Self::class_order`], keep all intermediate work in class order, and
/// pay for exactly one [`Self::public_coefficients`] epilogue at the end.
pub trait ClassOrderedCommutation<T, U> {
    /// The basis' class-contiguous ordering: created on first request
    /// (then cached on the series' table and shared), or prebuilt at
    /// series construction when a degree slice outgrows L1.
    fn class_order(&self) -> Arc<ClassOrder>;

    /// Rearranges public-order coefficients into class-contiguous order
    /// (the fold's operand conversion, once per fold).
    fn class_coefficients(&self, order: &ClassOrder, coefficients: &[U]) -> Vec<U>;

    /// The final epilogue: class-contiguous coefficients back to public
    /// basis order. Call once, after the last class-ordered kernel call.
    fn public_coefficients(&self, order: &ClassOrder, class_coefficients: &[U]) -> Vec<U>;

    /// Commutator entirely in class-contiguous order: `a`/`b` are
    /// class-ordered coefficient slices, the support lists class-indexed,
    /// and `result` receives the class-ordered sum directly — no scratch,
    /// no permutation. `order` must describe this series' basis.
    fn class_commutation(
        &self,
        order: &ClassOrder,
        a: &[U],
        a_nonzero: &[usize],
        b: &[U],
        b_nonzero: &[usize],
        result: &mut [U],
        cache: &mut GatingCache,
    );

    /// Batch form of [`Self::class_commutation`] for folds: jobs carry
    /// class-ordered operands, class-indexed support lists, and
    /// class-ordered result buffers.
    fn class_commutation_batch(
        &self,
        order: &ClassOrder,
        jobs: &mut [KernelJob<'_, U>],
        cache: &mut GatingCache,
    );

    /// Whole-fold form: `levels[l]` holds the jobs of dependency level `l`
    /// (operand slices may alias results of strictly earlier levels). The
    /// engine plans class-contiguous balanced packs for every level up front
    /// and runs the fold as one parallel pack-slot walk with level counters
    /// — no per-level dispatches.
    fn class_commutation_batch_fold(
        &self,
        order: &ClassOrder,
        levels: &mut [Vec<KernelJob<'_, U>>],
        cache: &mut GatingCache,
    ) {
        for level in levels.iter_mut() {
            if level.len() == 1 {
                self.class_commutation(
                    order,
                    // SAFETY: single-job serial path (exclusive).
                    unsafe { std::slice::from_raw_parts(level[0].a, level[0].a_len) },
                    level[0].a_nonzero,
                    unsafe { std::slice::from_raw_parts(level[0].b, level[0].b_len) },
                    level[0].b_nonzero,
                    unsafe { std::slice::from_raw_parts_mut(level[0].result, level[0].result_len) },
                    cache,
                );
            } else {
                self.class_commutation_batch(order, level, cache);
            }
        }
    }

    /// Collecting variant used by folds to (re)derive each result's
    /// scatter-target superset: the reported indices are class-indexed.
    #[allow(clippy::too_many_arguments)]
    fn class_commutation_with_nonzero_collecting(
        &self,
        order: &ClassOrder,
        a: &[U],
        a_nonzero: &[usize],
        b: &[U],
        b_nonzero: &[usize],
        result: &mut [U],
        dirty: &mut [u64],
        targets: &mut Vec<usize>,
    );
}

impl<T, U> ClassOrderedCommutation<T, U> for LieSeries<T, U>
where
    T: Clone + Ord + Generator + Hash,
    U: Clone
        + Default
        + One
        + Zero
        + Eq
        + MulAssign
        + Neg<Output = U>
        + Hash
        + Mul<Output = U>
        + AddAssign
        + Send
        + Sync
        + 'static,
{
    fn class_order(&self) -> Arc<ClassOrder> {
        self.feasible_decompositions.ensure_class_order()
    }

    fn class_coefficients(&self, order: &ClassOrder, coefficients: &[U]) -> Vec<U> {
        let inv = order.inv();
        let mut out = vec![U::default(); coefficients.len()];
        for (w, c) in coefficients.iter().enumerate() {
            out[inv[w] as usize] = c.clone();
        }
        out
    }

    fn public_coefficients(&self, order: &ClassOrder, class_coefficients: &[U]) -> Vec<U> {
        let inv = order.inv();
        let mut out = vec![U::default(); class_coefficients.len()];
        for (k, &src) in inv.iter().enumerate() {
            out[k] = class_coefficients[src as usize].clone();
        }
        out
    }

    fn class_commutation(
        &self,
        order: &ClassOrder,
        a: &[U],
        a_nonzero: &[usize],
        b: &[U],
        b_nonzero: &[usize],
        result: &mut [U],
        cache: &mut GatingCache,
    ) {
        let mut job = KernelJob {
            a: a.as_ptr(),
            a_len: a.len(),
            a_nonzero,
            b: b.as_ptr(),
            b_len: b.len(),
            b_nonzero,
            result: result.as_mut_ptr(),
            result_len: result.len(),
            a_shift: IDENTITY_SHIFTS.as_ptr(),
            b_shift: IDENTITY_SHIFTS.as_ptr(),
            r_shift: IDENTITY_SHIFTS.as_ptr(),
        };
        commutator_coefficients_class_batch_with_cache(
            self,
            order,
            std::slice::from_mut(&mut job),
            cache,
        );
    }

    fn class_commutation_batch(
        &self,
        order: &ClassOrder,
        jobs: &mut [KernelJob<'_, U>],
        cache: &mut GatingCache,
    ) {
        commutator_coefficients_class_batch_with_cache(self, order, jobs, cache);
    }

    fn class_commutation_with_nonzero_collecting(
        &self,
        order: &ClassOrder,
        a: &[U],
        a_nonzero: &[usize],
        b: &[U],
        b_nonzero: &[usize],
        result: &mut [U],
        dirty: &mut [u64],
        targets: &mut Vec<usize>,
    ) {
        LieSeries::commutator_coefficients_class_with_nonzero_collecting(
            self, order, a, a_nonzero, b, b_nonzero, result, dirty, targets,
        );
    }

    fn class_commutation_batch_fold(
        &self,
        order: &ClassOrder,
        levels: &mut [Vec<KernelJob<'_, U>>],
        cache: &mut GatingCache,
    ) {
        // The one-dispatch pack-slot walk is the default: its packs are
        // claimed dynamically, so a waiter never blocks on runnable fold
        // work and the walk stays live even when the pool is shared with
        // outer parallelism. Two cases route elsewhere: with a single
        // worker nothing can run in parallel, so the walk's planning (task
        // lists, packs, claim counters) would be pure overhead — the
        // per-level dispatch's serial branch sweeps the identical order
        // with no planning; and `FOLD_ENGINE=0` keeps the per-level
        // dispatch for A/B comparison.
        static PER_LEVEL_DISPATCH: OnceLock<bool> = OnceLock::new();
        let per_level = *PER_LEVEL_DISPATCH
            .get_or_init(|| std::env::var_os("FOLD_ENGINE").is_some_and(|v| v == *"0"));
        if per_level || rayon::current_num_threads() == 1 {
            for level in levels.iter_mut() {
                self.class_commutation_batch(order, level, cache);
            }
        } else {
            commutator_coefficients_class_fold_with_cache(self, order, levels, cache);
        }
    }
}

/// One stage of a class-ordered fold schedule: a kernel sweep level (packs
/// of `(job, unit)` tasks over one level's jobs) or a caller-defined block
/// stage (one atomic claim per block index, dispatched through a shared
/// closure — the gather/accumulate phases of a batched fold). A batch of
/// folds is a linear chain of stages; every slot waits at every stage
/// boundary, so a stage's reads observe exactly the writes of all earlier
/// stages (each completion is published with `Release`, observed with
/// `Acquire`).
pub struct ClassBatchStage<'a, U> {
    /// Sweep stages: `(job index within the stage, unit index in that
    /// job's gating)` tasks in sweep order. Block stages: unused.
    tasks: Arc<[(u32, u32)]>,
    /// Balanced pack ranges into `tasks` (sweep stages) or one block per
    /// pack (block stages).
    packs: Arc<[(usize, usize)]>,
    /// Sweep stages: per-job gating.
    gateways: Vec<KernelGating>,
    /// Sweep stages: the stage's jobs, indexed by the tasks' job ids.
    jobs: &'a [KernelJob<'a, U>],
    /// Block stages: `block(block_index, fold_index)` per claimed task.
    /// Blocks write disjoint ranges (the caller's contract); the stage
    /// counters order all cross-stage dependencies.
    block: Option<&'a (dyn Fn(usize, usize) + Send + Sync)>,
    /// Block stages: which fold of the batch this stage serves (the
    /// block task reads the fold's inputs through it). Sweep stages: 0.
    fold: usize,
    /// Sweep stages: run the 4-lane SoA cohort sweep instead of the scalar
    /// one. A cohort stage's jobs carry SoA-interleaved operand/result
    /// pointers (the four lanes' values of class position `p` live at
    /// `[p*4 .. p*4+4]`, lane `l` at `p*4 + l`), and every other plan
    /// input (tasks, packs, gateways, support lists) is lane-invariant —
    /// see `CommutatorDag::fold_batch_cohort` for how such a plan is
    /// built.
    cohort: bool,
}

/// The per-block task/pack shape shared by every block stage of one
/// batch: `blocks` disjoint-range tasks, one pack per task. Built once
/// per batch, then each block stage clones the `Arc`s (refcount bump —
/// no per-stage allocation).
pub struct BlockShape {
    tasks: Arc<[(u32, u32)]>,
    packs: Arc<[(usize, usize)]>,
}

impl BlockShape {
    /// The block shape over `blocks` disjoint-range tasks.
    pub fn new(blocks: usize) -> Self {
        Self {
            tasks: (0..blocks as u32).map(|i| (i, 0)).collect(),
            packs: (0..blocks).map(|i| (i, i + 1)).collect(),
        }
    }
}

impl<'a, U> ClassBatchStage<'a, U> {
    /// A block stage over `shape`'s disjoint-range tasks: pack `p` runs
    /// `task(p, fold)` exactly once.
    pub fn block_stage(
        shape: &BlockShape,
        task: &'a (dyn Fn(usize, usize) + Send + Sync),
        fold: usize,
    ) -> Self {
        Self {
            tasks: Arc::clone(&shape.tasks),
            packs: Arc::clone(&shape.packs),
            gateways: Vec::new(),
            jobs: &[],
            block: Some(task),
            fold,
            cohort: false,
        }
    }

    /// The cohort (SoA) sweep variant of a planned sweep stage: identical
    /// tasks, packs, gateways and job views (all lane-invariant — see the
    /// `cohort` field), dispatched to the 4-lane SoA kernel instead of the
    /// scalar sweep. The stage's jobs' operand/result pointers must be
    /// SoA-interleaved with lane stride 4 and the plan must have been
    /// built with `cohort = true` — see `plan_class_sweep_stages`.
    pub fn cohort_variant(stage: &Self) -> Self {
        assert!(
            stage.cohort,
            "cohort_variant requires a plan built with cohort = true              (SoA-interleaved buffers)"
        );
        Self {
            tasks: Arc::clone(&stage.tasks),
            packs: Arc::clone(&stage.packs),
            gateways: stage.gateways.clone(),
            jobs: stage.jobs,
            block: None,
            fold: 0,
            cohort: true,
        }
    }
}




/// A park/wake gate for one stage boundary. Waiters poll the stage
/// counter briefly (stage work is microsecond-scale; polling costs the
/// finisher nothing), and only a waiter whose gate stays closed parks —
/// registering itself so the finisher knows a kernel wake is warranted.
/// A finisher whose gate has no parked waiters pays a plain store.
/// (`done` remains the counter — it is what waiters check — the word only
/// carries the wake signal; `futex_wait` re-validates the word's value
/// atomically under the kernel lock, and the bounded park timeout
/// self-heals the benign race where a waiter registers just after the
/// finisher looked.)
struct FutexGate {
    word: AtomicU32,
    parked: AtomicUsize,
}

impl FutexGate {
    const fn new() -> Self {
        FutexGate {
            word: AtomicU32::new(0),
            parked: AtomicUsize::new(0),
        }
    }

    /// Sleeps while this gate's word is still `0`, bounded — a signal that
    /// raced past the wake check self-heals on timeout.
    fn park(&self) {
        self.parked.fetch_add(1, Ordering::AcqRel);
        loop {
            if self.word.load(Ordering::Acquire) != 0 {
                break;
            }
            #[cfg(target_os = "linux")]
            unsafe {
                let timeout = libc::timespec {
                    tv_sec: 0,
                    tv_nsec: 200_000, // 200 µs
                };
                let _ = libc::syscall(
                    libc::SYS_futex,
                    self.word.as_ptr(),
                    libc::FUTEX_WAIT | libc::FUTEX_PRIVATE_FLAG,
                    0u32,
                    &timeout,
                    0u32,
                    0u32,
                );
            }
            #[cfg(not(target_os = "linux"))]
            std::thread::sleep(std::time::Duration::from_micros(200));
        }
        self.parked.fetch_sub(1, Ordering::AcqRel);
    }

    /// Signals completion; wakes parked waiters if there are any.
    fn finish(&self) {
        self.word.store(1, Ordering::Release);
        if self.parked.load(Ordering::Acquire) > 0 {
            #[cfg(target_os = "linux")]
            unsafe {
                let _ = libc::syscall(
                    libc::SYS_futex,
                    self.word.as_ptr(),
                    libc::FUTEX_WAKE | libc::FUTEX_PRIVATE_FLAG,
                    i32::MAX,
                );
            }
        }
    }
}

/// Publishes one completed pack of `target`; the finisher that completes
/// the stage signals its gate. Returns true iff this call was the last
/// publisher.
#[inline]
fn finish_pack(done: &AtomicUsize, gate: &FutexGate, target: usize) -> bool {
    // AcqRel: the read half acquires the earlier finishers' releases (so
    // the last finisher's completion observes the completed stage), and
    // the write half releases this pack's writes to the next stage's readers.
    done.fetch_add(1, Ordering::AcqRel) + 1 == target && {
        gate.finish();
        true
    }
}

/// The batch walk's read-only view: the stage chain plus the sweep's
/// class-ordered tables. Operand slices ride inside the stages' jobs.
struct FoldWalk<'a, U> {
    stages: &'a [&'a ClassBatchStage<'a, U>],
    entries: &'a [Entry],
    /// Class position → Lyndon degree (the compact-buffer shifts index by
    /// degree; a contribution's entry positions are class positions).
    degrees: &'a [u8],
    coeffs: &'a [U],
    /// The class-space decomposition rel slice (rows are addressed by
    /// `(decomp_start, decomp_start+1)` spans into this slice; the
    /// coefficients share its indexing).
    decomp_rels: &'a [u32],
}

impl<'a, U: Clone + Neg<Output = U> + Mul<Output = U> + AddAssign + std::hash::Hash + 'static>
    FoldWalk<'a, U>
{
    /// Claims the next unclaimed pack of stage `s` off the stage's atomic
    /// cursor and runs it (a sweep pack or one block); `false` means the
    /// stage is fully claimed and its remaining work is in flight on other
    /// slots. Claims go to whoever asks, so a pack whose original runner is
    /// queued behind pool contention is run by any running slot — a waiter
    /// never waits on a runnable pack. (`Relaxed` suffices: claims only
    /// hand out disjoint pack indices; the produced data itself is
    /// published by the caller's `done` increment with `Release`.)
    #[inline]
    fn claim_and_run_pack(&self, s: usize, claim: &AtomicUsize) -> bool {
        let stage = self.stages[s];
        // Read-only exhaustion peek first: the gate waits drain this cursor
        // in a loop, and a *failed* claim must not RMW the line — 32 slots
        // spin-claiming an exhausted cursor would otherwise keep the cache
        // line exclusive-bouncing for the whole gate wait. Only a claim
        // that looks available touches the cursor.
        if claim.load(Ordering::Relaxed) >= stage.packs.len() {
            return false;
        }
        let pack = claim.fetch_add(1, Ordering::Relaxed);
        if pack >= stage.packs.len() {
            return false;
        }
        match stage.block {
            // Block stage: pack `p` runs block `p` (one task per pack).
            Some(task) => task(stage.tasks[pack].0 as usize, stage.fold),
            None if stage.cohort => self.cohort_sweep_pack_range(stage, pack),
            None => self.sweep_pack_range(stage, pack),
        }
        true
    }

    /// Sweeps pack `pack` of a sweep stage: single-phase per-unit
    /// direct-add. Packs are single-writer regions — a unit's word set
    /// partitions the write slots (ownership invariant) and packs never
    /// split a unit — so the per-word accumulation order is exactly the
    /// serial sweep's regardless of which slot claims which pack.
    #[inline]
    fn sweep_pack_range(&self, stage: &ClassBatchStage<'a, U>, pack: usize) {
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
            // SAFETY: within a stage the jobs' result buffers are pairwise
            // disjoint (caller contract), and within one job each unit's
            // word set is a single-writer region (two units never write the
            // same word, packs never split a unit), so concurrent packs of
            // the same job read-modify-write disjoint words through raw
            // pointers — no `&mut` is created, because two packs of one job
            // would otherwise hold aliasing `&mut` slices (UB under the
            // aliasing model even with disjoint ranges). Cross-stage
            // accesses are ordered by the stage counters. The buffer
            // outlives the walk (stable allocations, no resize during the
            // dispatch).
            let result: *mut U = job.result;
            let gating = &stage.gateways[ji];
            let a_shift = job.a_shift;
            let b_shift = job.b_shift;
            let r_shift = job.r_shift;
            for &(ji_task, ui) in &stage.tasks[t..e] {
                let _ = ji_task;
                let unit = &gating.units[ui as usize];
                let td = unit.td as usize;
                // The unit's degree-`td` slice starts at `rs` in the
                // working layout; the job's compact buffer stores the
                // target slice at `rs - r_shift[td]`.
                // SAFETY: the shift tables have >= max_degree+1 entries
                // and `td` is a degree of the table (see the plan).
                let base = unit.rs as usize - unsafe { *r_shift.add(td) } as usize;
                for tp in unit.ticket_start as usize..unit.ticket_end as usize {
                    let ticket = gating.tickets[tp];
                    let e_idx = (ticket & TICKET_INDEX_MASK) as usize;
                    let entry = self.entries[e_idx];
                    let p_active = ticket & TICKET_P_ACTIVE != 0;
                    let q_active = ticket & TICKET_Q_ACTIVE != 0;
                    let (i, j) = (entry.i as usize, entry.j as usize);
                    // Compact-layout operand shifts (see `KernelJob`): the
                    // contribution's operand degrees come from the
                    // class-position degree table.
                    let p = self.degrees[i] as usize;
                    let q = self.degrees[j] as usize;
                    // SAFETY: the shift tables have >= max_degree+1 entries
                    // and p/q are degrees of the table; i and j are class
                    // positions in the operands' supports (the tickets
                    // resolved presence). `raw_mul` skips the float
                    // wrappers' per-op NaN checks (raw-float fast path).
                    let term = unsafe {
                        if p_active {
                            let a_sh_p = *a_shift.add(p) as usize;
                            let b_sh_q = *b_shift.add(q) as usize;
                            let mut t =
                                raw_mul(&*job.a.add(i - a_sh_p), &*job.b.add(j - b_sh_q));
                            if q_active {
                                let a_sh_q = *a_shift.add(q) as usize;
                                let b_sh_p = *b_shift.add(p) as usize;
                                raw_add_assign(
                                    &mut t,
                                    &-raw_mul(&*job.a.add(j - a_sh_q), &*job.b.add(i - b_sh_p)),
                                );
                            }
                            t
                        } else {
                            let a_sh_q = *a_shift.add(q) as usize;
                            let b_sh_p = *b_shift.add(p) as usize;
                            -raw_mul(&*job.a.add(j - a_sh_q), &*job.b.add(i - b_sh_p))
                        }
                    };
                    let from = entry.decomp_start as usize;
                    let to = self.entries[e_idx + 1].decomp_start as usize;
                    for (k, &rel) in self.decomp_rels[from..to].iter().enumerate() {
                        // SAFETY: single-writer RMW (unit word-set
                        // ownership, packs never split a unit), in bounds
                        // (the gating addresses the job's compact result
                        // space; `rs + rel - r_shift[td]` is the compact
                        // slot of the row's target word).
                        unsafe {
                            raw_add_assign(
                                &mut *result.add(base + rel as usize),
                                &raw_mul(&self.coeffs[from + k], &term),
                            );
                        }
                    }
                }
            }
            t = e;
        }
    }

    /// The cohort (4-lane SoA) sweep for one pack — dispatched to the
    /// lane-typed kernel in [`cohort`] (the stage's `cohort` flag routes
    /// here; the engine that built the stage verified the cohort
    /// capability for `U`). The stage's jobs carry SoA-interleaved
    /// operand/result/term pointers; every other plan input (tasks, packs,
    /// gateways, support lists) is lane-invariant — see the `cohort`
    /// module's docs for the layout and the bit-exactness argument.
    #[inline]
    fn cohort_sweep_pack_range(&self, stage: &ClassBatchStage<'a, U>, pack: usize) {
        cohort::sweep_pack_range(
            self.entries,
            self.degrees,
            self.coeffs,
            self.decomp_rels,
            stage,
            pack,
        );
    }
}

/// Plans the sweep stages for the DAG's dependency levels: per-job gating
/// from the (fixed) support lists, then ONE single-phase per-unit sweep
/// stage per level (tasks = units, balanced pack cuts weighted by unit
/// pair counts). `cohort` builds the 4-lane SoA variant. The stages
/// borrow `levels`; a batched fold plans once and reuses the same stages
/// for every fold in the batch.
///
/// `want_scatter_sets` gates the exact per-job batch scatter sets: the
/// batch path sizes its compact buffers from them, but the per-fold path
/// discards them — computing them there is pure overhead (a `seen` bitvec
/// plus a full entry pass per job, O(jobs × (space + entries)) per level) —
/// so it passes `false` and receives an empty vec.
pub fn plan_class_sweep_stages<'a, T, U>(
    a_series: &LieSeries<T, U>,
    order: &ClassOrder,
    levels: &'a [Vec<KernelJob<'a, U>>],
    cache: &mut GatingCache,
    want_scatter_sets: bool,
    cohort: bool,
) -> (Vec<ClassBatchStage<'a, U>>, Vec<Vec<Vec<u32>>>)
where
    T: Clone + Ord + Generator + Hash,
    U: Clone
        + Default
        + One
        + Zero
        + Eq
        + MulAssign
        + Neg<Output = U>
        + Hash
        + Mul<Output = U>
        + AddAssign
        + Send
        + Sync
        + 'static,
{
    let threads = rayon::current_num_threads().max(1);
    let mut stages = Vec::with_capacity(levels.len());
    // The jobs' TRUE batch scatter sets: the class positions each node's
    // sweep writes under its (fixed) operand supports. The recorded node
    // lists are only a union-level fixed point (the batch eligibility
    // checks the union, not each list) — a level-0 job's list, recorded
    // under an earlier, smaller accumulator support, can miss positions
    // the batch's first fold scatters. Compact buffers must therefore be
    // sized from this exact set. Under the per-unit division the set is
    // exact and free: the gating's `unit_words` list IS the write set
    // (ascending — units in degree order, word sets ascending per unit).
    let mut scatter_sets: Vec<Vec<Vec<u32>>> =
        Vec::with_capacity(if want_scatter_sets { levels.len() } else { 0 });
    // Pass 1: gate every job up front (the support lists are fixed for the
    // whole fold/batch, so no stage's gating depends on kernel results)
    // and total the planned sweep work — the fold-level slot budget that
    // every stage's pack cut below shares. A stage pack exists to give a
    // walk slot work during that stage, so packs must never outnumber the
    // slots the policy will actually run: with the budget derived from the
    // whole fold's work, a tiny fold (12x2: ~66 entries) plans ONE pack per
    // stage and the serial walk drains a single claim instead of nine.
    let level_gateways: Vec<Vec<KernelGating>> = levels
        .iter()
        .map(|level| {
            level
                .iter()
                .map(|j| {
                    LieSeries::<T, U>::kernel_prologue_cached_class(
                        a_series,
                        j.a_nonzero,
                        j.b_nonzero,
                        order,
                        cache,
                    )
                })
                .collect()
        })
        .collect();
    // Planned sweep work per fold: one term per active entry PLUS one
    // multiply-add per contribution pair (the single-phase unit sweep's
    // two work units — the same magnitudes the two-phase planner funded).
    let sweep_entries: usize = level_gateways
        .iter()
        .map(|gateways| {
            gateways
                .iter()
                .map(|g| g.total_pairs + g.tickets.len())
                .sum::<usize>()
        })
        .sum();
    let slots = work_adaptive_slots(threads, usize::MAX, sweep_entries);
    for (level, gateways) in levels.iter().zip(level_gateways) {
        if want_scatter_sets {
            // Exact scatter sets for this level's jobs: the active word
            // positions themselves (each is written exactly once by the
            // sweep, in ascending order).
            let level_sets: Vec<Vec<u32>> = gateways
                .iter()
                .map(|gateway| gateway.unit_words.to_vec())
                .collect();
            scatter_sets.push(level_sets);
        }
        // Flatten (job, unit) tasks: the per-unit division's atomic
        // scheduling units. Units are emitted per job in table (degree)
        // order and packs cut only between units — race-free (unit word
        // sets partition the write slots) and bit-exact (a word's adds
        // stay inside its one unit, in table-entry order).
        let total: usize = gateways.iter().map(|g| g.total_pairs).sum();
        let unit_count: usize = gateways.iter().map(|g| g.units.len()).sum();
        let mut tasks: Vec<(u32, u32)> = Vec::with_capacity(unit_count);
        for (ji, gating) in gateways.iter().enumerate() {
            for ui in 0..gating.units.len() {
                tasks.push((ji as u32, ui as u32));
            }
        }
        // Balanced prefix cuts at unit boundaries: pack `p` of `p_count`
        // takes tasks while `cur < total * (p + 1) / p_count`. Packs are
        // capped by the task count, not the job count: a narrow upper
        // level's few big jobs must still spread across every slot (unit
        // word sets are conflict-free write regions, and each word's
        // accumulation sequence stays whole inside whichever pack owns its
        // unit).
        //
        // The balance weight is the unit's pair count (its WORK proxy: one
        // term computation per active entry plus one multiply-add per
        // pair).
        // Smallest work worth splitting a stage's sweep across an extra
        // pack: below ~8 contributions per pack the per-pack claim/gate
        // cost dominates the pack's compute.
        const MIN_PACK_WORK: usize = 8;
        let p_count = slots
            .min(unit_count.max(1))
            // Work-scaled cap: a stage too small to give every slot a
            // pack of at least `MIN_PACK_WORK` contributions must not
            // spawn micro-packs — the claim/gate churn outweighs the
            // parallelism on narrow grids (e.g. 12x2, whose whole fold is
            // ~66 entries per stage). Balanced bounds below keep the
            // resulting packs even. The leading `slots` cap is the
            // fold-level work budget (see pass 1): packs beyond the
            // policy's slot count would only be drained serially.
            .min(total.div_ceil(MIN_PACK_WORK))
            .max(1);
        let mut packs: Vec<(usize, usize)> = Vec::with_capacity(p_count);
        if p_count > 1 && !tasks.is_empty() {
            let mut start = 0usize;
            let mut cur = 0usize;
            for (ti, &(ji, ui)) in tasks.iter().enumerate() {
                let unit = &gateways[ji as usize].units[ui as usize];
                cur += unit.pairs as usize;
                let next_bound = total * (packs.len() + 1) / p_count;
                if cur >= next_bound && ti + 1 < tasks.len() {
                    packs.push((start, ti + 1));
                    start = ti + 1;
                }
            }
            packs.push((start, tasks.len()));
        } else if !tasks.is_empty() {
            packs.push((0, tasks.len()));
        }
        stages.push(ClassBatchStage {
            tasks: tasks.into(),
            packs: packs.into(),
            gateways,
            jobs: level.as_slice(),
            block: None,
            fold: 0,
            cohort,
        });
    }
    (stages, scatter_sets)
}

/// The slot policy's work quantum, in planned active-entry tickets per
/// fold unit per slot: one more coordinated slot must be funded by this
/// much real sweep work. Calibration: a swept entry costs ~40 cycles
/// (measured marginal on the class kernel) and a coordinated slot pays
/// the per-stage claim/gate protocol (~2-5 µs across a stage boundary),
/// so 3750 entries × 40 cy ≈ 150 Kcy ≈ 50 µs of sweep work per slot at
/// ~3 GHz — the break-even where an extra slot stops paying for the
/// barriers it joins. Tuned between the two candidates (1500 ≈ 20 µs and
/// 3750 ≈ 50 µs equivalents): both fix the tiny regimes, 3750 also wins
/// the one regime they disagree on (3x8 at a 32t pool: e2e 161 ms vs
/// 171 ms). Measured regimes: 12x2 (66 entries/fold) and 8x3 (~700/fold)
/// walk serially at every pool size (their e2e used to run 6.5x slower at
/// 16t than 1t — barrier machinery dwarfing the work); 3x8 (~74K/fold)
/// lands near its measured 8-16t sweet spot (19 slots at a 32t pool);
/// 2x12 (~247K/fold) still fills the pool.
const SLOT_WORK_QUANTUM: usize = 3750;

/// SLOT_POLICY_DEBUG=1: one stderr line per walk dispatch with the
/// policy's inputs and decision (tuning/ops diagnostics; off by default).
fn slot_policy_debug() -> bool {
    static DEBUG: OnceLock<bool> = OnceLock::new();
    *DEBUG.get_or_init(|| std::env::var_os("SLOT_POLICY_DEBUG").is_some())
}

/// The stage chain's planned SWEEP work in swept-entry units: sweep
/// stages contribute their gateways' precomputed active-entry tickets (the
/// exact per-entry work unit the sweeps iterate — presence-resolved, so
/// canceled entries don't inflate the estimate); block stages contribute
/// nothing (their element count depends on the caller's block cut, which
/// is itself chosen from this estimate). Pure plan-time data: no kernel
/// results are needed.
pub fn planned_sweep_entries<U>(stages: &[&ClassBatchStage<'_, U>]) -> usize {
    stages
        .iter()
        .map(|s| match s.block {
            Some(_) => 0,
            None => s
                .gateways
                .iter()
                .map(|g| g.total_pairs + g.tickets.len())
                .sum::<usize>(),
        })
        .sum()
}

/// The stage chain's planned work in swept-entry units: sweep stages
/// contribute their gateways' precomputed active-entry tickets, block
/// stages their element count (one block = one disjoint class-position
/// range per pack). Pure plan-time data: no kernel results are needed.
fn planned_stage_entries<U>(stages: &[&ClassBatchStage<'_, U>]) -> usize {
    stages
        .iter()
        .map(|s| match s.block {
            Some(_) => s.packs.len(),
            None => s
                .gateways
                .iter()
                .map(|g| g.total_pairs + g.tickets.len())
                .sum::<usize>(),
        })
        .sum()
}

/// Work-adaptive slot count for the stage-chain walk: the walk's parallel
/// width is chosen from the PLANNED work, not just the pool size. Every
/// slot joins every stage boundary's counter/gate protocol, so a slot is
/// worth spawning only when its share of the per-fold-unit work covers the
/// coordination cost: `per_unit_work / QUANTUM` slots, capped by the pool
/// and by the widest stage's pack count (a slot beyond every stage's packs
/// only adds barrier arrivals, never runs work). Callers that cut their
/// own stage structure (e.g. a batch's block count) should derive it from
/// the same policy so packs never outnumber useful slots.
pub fn work_adaptive_slots(threads: usize, max_packs: usize, per_unit_work: usize) -> usize {
    threads
        .min(max_packs)
        .min((per_unit_work / SLOT_WORK_QUANTUM).max(1))
        .max(1)
}

/// Runs a linear chain of fold stages as ONE parallel dispatch (a serial
/// loop on the calling thread whenever the planned work cannot fund more
/// than one coordinated slot). Every slot walks the stages in order:
/// before waiting on a stage's counter it drains that stage's claim cursor
/// — a queued stage's packs are run by whoever is working, never waited
/// for — so the walk stays live under any pool contention. The per-stage
/// counters (Release on completion, Acquire at the boundary) are the only
/// synchronization: a stage's reads observe exactly the writes of all
/// earlier stages, and the per-word accumulation order is exactly the
/// serial schedule's (packs never split a unit; blocks never split a
/// range).
///
/// `fold_units` is how many per-fold repetitions the stage chain contains
/// (1 for the per-fold path; a batch of `rhss.len()` folds repeats its
/// gather/sweep/accumulate sub-chain once per displacement) — the slot
/// policy normalizes the planned work by it, because the barrier cost is
/// paid per stage, i.e. per fold unit, not once per dispatch.
pub fn run_class_batch<T, U>(
    a_series: &LieSeries<T, U>,
    order: &ClassOrder,
    stages: &[&ClassBatchStage<'_, U>],
    fold_units: usize,
) where
    T: Clone + Ord + Generator + Hash,
    U: Clone
        + Neg<Output = U>
        + Mul<Output = U>
        + AddAssign
        + std::hash::Hash
        + Send
        + Sync
        + 'static,
{
    run_class_batch_with_work(
        a_series,
        order,
        stages,
        fold_units,
        planned_stage_entries(stages),
    );
}

/// [`run_class_batch`] with the caller's planned-work estimate: a cohort
/// chain's sweep stages carry the 4-lane work that the plain stage-based
/// estimate (entry tickets, lane-agnostic plan data) cannot see, so the
/// cohort engine supplies `planned_work = scalar_tickets × lanes` and the
/// slot policy funds the walk from the true per-unit work.
pub fn run_class_batch_with_work<T, U>(
    a_series: &LieSeries<T, U>,
    order: &ClassOrder,
    stages: &[&ClassBatchStage<'_, U>],
    fold_units: usize,
    planned_work: usize,
) where
    T: Clone + Ord + Generator + Hash,
    U: Clone
        + Neg<Output = U>
        + Mul<Output = U>
        + AddAssign
        + std::hash::Hash
        + Send
        + Sync
        + 'static,
{
    // The walk is fully internal for sweep stages: relabeled entries gate
    // the presence tests and index the class-ordered operands. Results are
    // written per job through the job's own buffer (see `sweep_pack_range`).
    let walk = FoldWalk {
        stages,
        entries: order.entries_cls(),
        degrees: order.degree_cls(),
        coeffs: a_series.feasible_decompositions.decomp_coeffs(),
        decomp_rels: order.decomp_cls(),
    };

    // Pack-slot walk with dynamic claims: `slots` tasks walk every stage,
    // and each stage's packs are handed out by an atomic cursor at run
    // time. A slot therefore never waits on a *queued* runner: before it
    // blocks on a stage counter it drains that stage's claim cursor itself,
    // so every unclaimed pack is either claimed-and-in-flight (a running
    // task publishes it via the counter) or claimable by the waiting slot.
    // That keeps the walk live even when the pool is shared with outer
    // parallelism and slot tasks queue behind other work — the fixed
    // slot→pack assignment this replaces starved exactly there.
    let stage_pack_counts: Vec<usize> = stages.iter().map(|s| s.packs.len()).collect();
    // One padded cell per counter: adjacent stages' counters must not share
    // a cache line — 32 slots RMW-ing two hot counters per gate otherwise
    // drag every *other* stage's counters along with them.
    #[repr(align(128))]
    struct Pad(AtomicUsize);
    let done: Vec<Pad> = stage_pack_counts
        .iter()
        .map(|_| Pad(AtomicUsize::new(0)))
        .collect();
    let claims: Vec<Pad> = stage_pack_counts
        .iter()
        .map(|_| Pad(AtomicUsize::new(0)))
        .collect();
    let gates: Vec<FutexGate> = stage_pack_counts.iter().map(|_| FutexGate::new()).collect();
    let threads = rayon::current_num_threads().max(1);
    let max_packs = stage_pack_counts.iter().copied().max().unwrap_or(1);
    let planned_entries = planned_work;
    let per_unit_work = planned_entries / fold_units.max(1);
    let slots = work_adaptive_slots(threads, max_packs, per_unit_work);
    if slot_policy_debug() {
        eprintln!(
            "slot_policy: fold_units={fold_units} planned_entries={planned_entries} \
             per_unit={per_unit_work} threads={threads} max_packs={max_packs} slots={slots}"
        );
    }
    let walk_for_slot = |_slot: usize| {
        // Hoisted per-stage sync references: the drain/spin loops below run
        // per claim attempt / per spin iteration, and re-indexing the
        // counter Vecs there cost ~7% of the wide-pool profile in bounds
        // checks + deref chains alone (flamegraph, 4x8/32t). The Vecs are
        // never resized during the walk, so shared references hoisted once
        // per slot are equivalent.
        let stage_done: Vec<&AtomicUsize> = done.iter().map(|p| &p.0).collect();
        let stage_claims: Vec<&AtomicUsize> = claims.iter().map(|p| &p.0).collect();
        let stage_gates: Vec<&FutexGate> = gates.iter().collect();
        for s in 0..stages.len() {
            if s > 0 {
                // EVERY slot waits at every stage boundary — including
                // slots that claimed nothing at the previous stage — so no
                // task of stage `s` can read a buffer written by stage
                // `s - 1` before all of stage `s - 1`'s packs have released
                // (Release in the run callers, Acquire here). While
                // waiting, keep draining the claim cursor: if the slots
                // that claimed the stage's remaining packs are queued
                // behind pool contention, this running slot runs those
                // packs itself instead of blocking on them. Spin briefly,
                // then yield to co-scheduled rayon work, then sleep — at
                // this point every remaining pack is in flight, so the
                // wait is bounded by real work.
                let need = stage_pack_counts[s - 1];
                let (prev_done, prev_gate, prev_claim) =
                    (stage_done[s - 1], stage_gates[s - 1], stage_claims[s - 1]);
                // Drain once on arrival: after this the cursor is exhausted
                // — every pack is claimed-and-in-flight — so parking is
                // safe (no pack waits on this slot). The cursor never
                // grows back, so no re-drain is ever needed.
                while walk.claim_and_run_pack(s - 1, prev_claim) {
                    finish_pack(prev_done, prev_gate, need);
                }
                // Poll briefly (most gates open within microseconds —
                // polling costs the finisher nothing), then YIELD instead
                // of parking: a parked slot BLOCKS its rayon worker, and
                // the tournament runs many concurrent dispatches on one
                // pool (e2e: 32 leaf-group tasks × slot tasks on 32
                // workers) — queued sibling slot tasks sit behind the
                // parked holder, which measured as 40% slot utilization
                // at 32 threads (task-clock vs wall, 4x8 concat). A
                // yielding slot re-queues its task so the worker runs the
                // other dispatches' pack work; the drain-on-arrival rule
                // guarantees every pack is claimed-and-in-flight, so the
                // gate always opens from real work. Parking stays as the
                // long-wait fallback (frees the core when a stage is
                // genuinely slow; the finisher then pays the wake).
                let mut spins = 0usize;
                let mut yields = 0usize;
                while prev_done.load(Ordering::Acquire) < need {
                    if spins < 1 << 12 {
                        std::hint::spin_loop();
                        spins += 1;
                    } else if yields < 1 << 10 {
                        yields += 1;
                        rayon::yield_now();
                    } else {
                        prev_gate.park();
                    }
                }
            }
            // Run this stage's packs off the claim cursor and publish each
            // completion: `done[s]` reaches exactly `packs[s]` when — and
            // only when — every pack has been run by whichever slot
            // claimed it. (A slot that finds the cursor drained reports
            // nothing, or the stage would appear complete before its work
            // is done.)
            let (this_done, this_gate, this_claim) =
                (stage_done[s], stage_gates[s], stage_claims[s]);
            let this_need = stage_pack_counts[s];
            while walk.claim_and_run_pack(s, this_claim) {
                finish_pack(this_done, this_gate, this_need);
            }
        }
    };

    // The serial walk is entered by the POLICY's slot count, not the pool
    // size: `slots == 1` runs the identical walk on the calling thread —
    // the claims/counters/gates are self-contained per slot, so the
    // single-slot semantics are exactly the single-worker pool's long-
    // standing serial path, minus the par_iter dispatch cost.
    if slots <= 1 {
        walk_for_slot(0);
    } else {
        use rayon::prelude::*;
        (0..slots).into_par_iter().for_each(walk_for_slot);
    }
}

/// Multi-level batch kernel in the class-contiguous working mode: the whole
/// fold runs as ONE parallel dispatch.
///
/// `levels[l]` holds the jobs of dependency level `l`; every job's operands
/// are class-ordered slices of buffers produced by strictly earlier levels
/// (or the fold's inputs), its support lists are class-indexed, and its
/// result buffer receives the class-ordered sum directly. Jobs of one level
/// must have pairwise disjoint result buffers.
///
/// The call plans every level up front (gating from the fixed support
/// lists, then balanced class-contiguous packs) and runs one parallel walk
/// of sweep stages. Each stage's packs are claimed dynamically off an
/// atomic cursor, every completion is published with `Release`, and every
/// slot waits (`Acquire`) at every stage boundary — so cross-level operand
/// reads are ordered by the stage counters, the per-word accumulation
/// order is exactly the serial sweep's (units are never split across
/// packs), and a slot never blocks on runnable fold work: before waiting
/// it drains the level's claim cursor itself.
pub fn commutator_coefficients_class_fold_with_cache<T, U>(
    a_series: &LieSeries<T, U>,
    order: &ClassOrder,
    levels: &mut [Vec<KernelJob<'_, U>>],
    cache: &mut GatingCache,
) where
    T: Clone + Ord + Generator + Hash,
    U: Clone
        + Default
        + One
        + Zero
        + Eq
        + MulAssign
        + Neg<Output = U>
        + Hash
        + Mul<Output = U>
        + AddAssign
        + Send
        + Sync
        + 'static,
{
    debug_assert_eq!(
        order.inv().len(),
        a_series.basis.len(),
        "class ordering does not describe this series' basis"
    );

    // The per-fold path discards the scatter sets (only the batch path
    // sizes compact buffers from them): skip their exact-set computation —
    // the gating and pack cuts the sweep actually reads are unchanged.
    let (stages, _scatter_sets) =
        plan_class_sweep_stages(a_series, order, levels, cache, false, false);
    let stage_refs: Vec<&ClassBatchStage<U>> = stages.iter().collect();
    // One fold unit: the stage chain is a single fold's sweep stages, so
    // the slot policy sees the fold's own planned work.
    run_class_batch(a_series, order, &stage_refs, 1);
}

/// Receiver for the absolute basis indices a kernel call scatters onto.
/// The collecting implementation deduplicates through a dirty bitset.
pub(crate) trait ScatterSink {
    #[inline(always)]
    fn scatter(&mut self, _index: usize) {}
}

pub(crate) struct NoSink;

impl ScatterSink for NoSink {}

pub(crate) struct CollectSink<'a> {
    dirty: &'a mut [u64],
    out: &'a mut Vec<usize>,
    cutoff: usize,
}

impl ScatterSink for CollectSink<'_> {
    #[inline(always)]
    fn scatter(&mut self, index: usize) {
        if index >= self.cutoff {
            return;
        }
        let word = &mut self.dirty[index / 64];
        let bit = 1u64 << (index % 64);
        if *word & bit == 0 {
            *word |= bit;
            self.out.push(index);
        }
    }
}

impl<
    T: Clone + Ord + Generator + Hash + Eq,
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
> LieSeries<T, U>
{
    /// Indices of non-zero coefficients that the commutation kernel will
    /// actually use: non-zero values on basis elements below `max_degree`.
    pub fn nonzero_coefficient_indices(&self, coefficients: &[U]) -> Vec<usize> {
        let cutoff = self
            .feasible_decompositions
            .degree_start(self.max_degree)
            .min(coefficients.len());
        coefficients[..cutoff]
            .iter()
            .enumerate()
            .filter(|(_, c)| !c.is_zero())
            .map(|x| x.0)
            .collect()
    }

    /// Total number of stored feasible decomposition pairs of this series'
    /// table — the static per-fold sweep-work upper bound (a steady fold's
    /// gating resolves most of the table's entries once operand supports
    /// are dense).
    pub fn feasible_table_len(&self) -> usize {
        self.feasible_decompositions.len()
    }

    /// Allocless [`Self::nonzero_coefficient_indices`] equality test:
    /// `true` iff `coefficients`' kernel-visible nonzero support equals
    /// `support` (both understood over `[0, cutoff)`, `support` sorted
    /// ascending — exactly what the reference call returns). Walks the
    /// coefficients once, matching nonzero positions against `support`
    /// in order; early-exits on the first mismatch. Same result as
    /// `self.nonzero_coefficient_indices(coefficients) == support`
    /// without the per-candidate `Vec` allocation.
    pub fn has_support(&self, coefficients: &[U], support: &[usize]) -> bool {
        let cutoff = self
            .feasible_decompositions
            .degree_start(self.max_degree)
            .min(coefficients.len());
        let mut remaining = support.iter().copied();
        let mut expected = remaining.next();
        for (i, c) in coefficients[..cutoff].iter().enumerate() {
            if !c.is_zero() {
                match expected {
                    Some(j) if j == i => expected = remaining.next(),
                    _ => return false,
                }
            }
        }
        expected.is_none()
    }

    #[cfg_attr(
        feature = "tracing",
        tracing::instrument(
            level = "debug",
            name = "lie_series_commutator_coefficients",
            skip_all,
            fields(
                basis_len = a_series.basis.len(),
                max_degree = a_series.max_degree
            )
        )
    )]
    pub fn commutator_coefficients(
        a_series: &LieSeries<T, U>,
        a_coefficients: &[U],
        b_coefficients: &[U],
        result_coefficients: &mut [U],
    ) where
        U: Send + Sync,
    {
        let a_nonzero = a_series.nonzero_coefficient_indices(a_coefficients);
        let b_nonzero = a_series.nonzero_coefficient_indices(b_coefficients);

        LieSeries::commutator_coefficients_with_nonzero(
            a_series,
            a_coefficients,
            &a_nonzero,
            b_coefficients,
            &b_nonzero,
            result_coefficients,
        );
    }

    pub fn commutator_coefficients_with_nonzero_collecting(
        a_series: &LieSeries<T, U>,
        a_coefficients: &[U],
        a_nonzero: &[usize],
        b_coefficients: &[U],
        b_nonzero: &[usize],
        result_coefficients: &mut [U],
        dirty: &mut [u64],
        targets: &mut Vec<usize>,
    ) {
        for word in dirty.iter_mut() {
            *word = 0;
        }
        targets.clear();
        let cutoff = a_series
            .feasible_decompositions
            .degree_start(a_series.max_degree);
        let gating = Self::kernel_prologue(a_series, a_nonzero, b_nonzero);
        let table = &a_series.feasible_decompositions;
        // The gating's word classes carry public positions: one single-writer
        // store per word, straight into the public result buffer (the old
        // class-scratch + permutation detour densified a scatter of repeated
        // `+=`s — one store per word needs no densifying). The sink reports
        // public basis indices directly.
        Self::sweep_words_serial(
            a_series,
            a_coefficients,
            b_coefficients,
            &gating,
            result_coefficients,
            table.entries(),
            table.decomp_indices_rel(),
            &mut CollectSink {
                dirty,
                out: targets,
                cutoff,
            },
        );
    }

    /// Class-contiguous variant of
    /// [`Self::commutator_coefficients_with_nonzero_collecting`]: operands
    /// are class-ordered, support lists class-indexed, and the sink reports
    /// class-indexed positions. No scratch, no permutation.
    #[allow(clippy::too_many_arguments)]
    /// Handle to the compiled feasible-decomposition table (basis-derived
    /// integer data — no `T`). `#[doc(hidden)]` kernel-plumbing accessor
    /// for the DAG crate's level-parallel collecting rebuild, which must
    /// capture only Send+Sync, T-free state (see `class_collect_kernel`).
    #[doc(hidden)]
    pub fn feasible_decompositions_handle(&self) -> &Arc<FeasibleDecompositions<U>> {
        &self.feasible_decompositions
    }

    pub(crate) fn commutator_coefficients_class_with_nonzero_collecting(
        a_series: &LieSeries<T, U>,
        order: &ClassOrder,
        a_coefficients: &[U],
        a_nonzero_cls: &[usize],
        b_coefficients: &[U],
        b_nonzero_cls: &[usize],
        result_coefficients: &mut [U],
        dirty: &mut [u64],
        targets: &mut Vec<usize>,
    ) {
        Self::class_collect_kernel(
            &a_series.feasible_decompositions,
            a_series.basis.len(),
            a_series.max_degree,
            order,
            a_coefficients,
            a_nonzero_cls,
            b_coefficients,
            b_nonzero_cls,
            result_coefficients,
            dirty,
            targets,
        );
    }

    /// T-free core of
    /// [`Self::commutator_coefficients_class_with_nonzero_collecting`]:
    /// everything the class collect kernel touches is basis-derived
    /// integer data (the feasible table, the class order) plus U values —
    /// the letter type never enters the kernel. This is what lets the
    /// collecting rebuild run level-parallel without shipping `T` across
    /// threads (see `CommutatorDag::ensure_lists_steady`).
    #[doc(hidden)]
    pub fn class_collect_kernel(
        table: &FeasibleDecompositions<U>,
        basis_len: usize,
        max_degree: usize,
        order: &ClassOrder,
        a_coefficients: &[U],
        a_nonzero_cls: &[usize],
        b_coefficients: &[U],
        b_nonzero_cls: &[usize],
        result_coefficients: &mut [U],
        dirty: &mut [u64],
        targets: &mut Vec<usize>,
    ) {
        for word in dirty.iter_mut() {
            *word = 0;
        }
        targets.clear();
        let cutoff = table.degree_start(max_degree);
        let gating = Self::kernel_prologue_cached_class_core(
            table,
            basis_len,
            a_nonzero_cls,
            b_nonzero_cls,
            order,
            &mut GatingCache::default(),
        );
        Self::sweep_words_serial_core(
            table.decomp_coeffs(),
            a_coefficients,
            b_coefficients,
            &gating,
            result_coefficients,
            order.entries_cls(),
            order.decomp_cls(),
            &mut CollectSink {
                dirty,
                out: targets,
                cutoff,
            },
        );
    }

    pub fn commutator_coefficients_with_nonzero(
        a_series: &LieSeries<T, U>,
        a_coefficients: &[U],
        a_nonzero: &[usize],
        b_coefficients: &[U],
        b_nonzero: &[usize],
        result_coefficients: &mut [U],
    ) where
        U: Send + Sync,
    {
        let mut job = KernelJob {
            a: a_coefficients.as_ptr(),
            a_len: a_coefficients.len(),
            a_nonzero,
            b: b_coefficients.as_ptr(),
            b_len: b_coefficients.len(),
            b_nonzero,
            result: result_coefficients.as_mut_ptr(),
            result_len: result_coefficients.len(),
            a_shift: IDENTITY_SHIFTS.as_ptr(),
            b_shift: IDENTITY_SHIFTS.as_ptr(),
            r_shift: IDENTITY_SHIFTS.as_ptr(),
        };
        commutator_coefficients_batch(a_series, std::slice::from_mut(&mut job));
        // `job` moved into the slice; its fields are all references/raw
        // pointers — nothing to drop.
    }

    fn kernel_prologue(
        a_series: &LieSeries<T, U>,
        a_nonzero: &[usize],
        b_nonzero: &[usize],
    ) -> KernelGating {
        let mut cache = GatingCache::default();
        Self::kernel_prologue_cached(a_series, a_nonzero, b_nonzero, &mut cache)
    }
    fn kernel_prologue_cached(
        a_series: &LieSeries<T, U>,
        a_nonzero: &[usize],
        b_nonzero: &[usize],
        cache: &mut GatingCache,
    ) -> KernelGating {
        let table = &a_series.feasible_decompositions;

        // Full-support fast path: the nonzero lists must be EXACTLY the
        // kernel-visible prefix `0..cutoff` — every index the sweep's
        // presence tests can ever reference (entries only pair positions
        // strictly below `max_degree`). Length alone is NOT sufficient:
        // node support lists recorded by a batch fold legitimately contain
        // degree-`max_degree` positions, and a same-length list that
        // includes them would skip the presence tests and read compact
        // slots the operand layout does not cover (wrong values — see the
        // pooled-dag leaf-reuse regression test).
        let cutoff = table.degree_start(table.max_degree());
        let full_prefix = |s: &[usize]| s.len() == cutoff && s.iter().all(|&p| p < cutoff);
        if full_prefix(a_nonzero) && full_prefix(b_nonzero) {
            let (units, tickets, unit_words) = table.full_support_gating_public();
            return KernelGating {
                total_pairs: units.iter().map(|u| u.pairs as usize).sum(),
                units,
                tickets,
                unit_words,
            };
        }

        // Memoized gating keyed by the exact support lists (see
        // `GatingCache`): on a hit neither the bitsets, the per-entry
        // presence tests, nor the per-unit word-set collection run at all.
        let key = (
            Self::support_fingerprint(a_nonzero),
            Self::support_fingerprint(b_nonzero),
        );
        if let Some((units, tickets, unit_words, total)) = cache.get(key) {
            return KernelGating {
                units: units.clone(),
                tickets: tickets.clone(),
                unit_words: unit_words.clone(),
                total_pairs: *total,
            };
        }
        // Fresh gating: presence bitsets drive the per-entry ticket
        // flags, degree-support masks the run-level gating. A unit's
        // entries are grouped into contiguous p-runs (entries are sorted
        // `(p, q, i, j)` within the unit and `q = target - p` is forced),
        // and a run is kept only when its own `(p, target - p)` degree
        // pair is supported — unit-level gating would drag every other p's
        // entries through tests that always fail. Surviving entries carry
        // pre-packed orientation flags, and the run walk simultaneously
        // collects each active entry's row targets into the owning unit's
        // word set, so no sweep re-derives anything.
        let words = a_series.basis.len().div_ceil(64);
        let mut presence = vec![0u64; 2 * words];
        let (a_present, b_present) = presence.split_at_mut(words);
        let mut a_deg = [0u64; 2];
        let mut b_deg = [0u64; 2];
        for &i in a_nonzero {
            a_present[i / 64] |= 1u64 << (i % 64);
            let d = table.degree_of(i);
            debug_assert!(d < 128, "degree masks cover degrees 0..127");
            a_deg[d / 64] |= 1u64 << (d % 64);
        }
        for &j in b_nonzero {
            b_present[j / 64] |= 1u64 << (j % 64);
            let d = table.degree_of(j);
            debug_assert!(d < 128, "degree masks cover degrees 0..127");
            b_deg[d / 64] |= 1u64 << (d % 64);
        }
        let value = Self::build_gating(
            table,
            table.entries(),
            table.entries(),
            table.decomp_indices_rel(),
            a_present,
            b_present,
            a_deg,
            b_deg,
        );
        let KernelGating {
            ref units,
            ref tickets,
            ref unit_words,
            total_pairs,
        } = value;
        cache.insert(
            key,
            (
                Arc::clone(units),
                Arc::clone(tickets),
                Arc::clone(unit_words),
                total_pairs,
            ),
        );
        value
    }

    /// 128-bit fingerprint of a support list, the gating cache's key half.
    /// Content-addressed: an entry's baked per-entry presence results are
    /// only valid for the exact supports they were resolved from, so two
    /// calls may share a cache entry only when their support lists are
    /// identical. Length is mixed in and both streams are avalanche-
    /// finished (murmur3 fmix64), so lists of different lengths or content
    /// collide with negligible probability.
    /// The gating cache's key fingerprint for one support list (murmur3-
    /// finish FNV pairs). Length is mixed in and both streams are
    /// avalanche-finished, so distinct lists collide with negligible
    /// probability.
    #[inline]
    fn support_fingerprint(list: &[usize]) -> [u64; 2] {
        let mut h1 = 0xcbf2_9ce4_8422_2325u64 ^ (list.len() as u64).rotate_left(32);
        let mut h2 = 0x9ae1_6a3b_2f90_404fu64 ^ (list.len() as u64);
        for &x in list {
            let v = x as u64;
            h1 = (h1 ^ v).wrapping_mul(0x1000_0000_01b3);
            h2 = (h2 ^ v).wrapping_mul(0x94d0_49bb_1331_11eb);
        }
        h1 ^= h1 >> 33;
        h1 = h1.wrapping_mul(0xff51_afd7_ed55_8ccd);
        h1 ^= h1 >> 33;
        h2 ^= h2 >> 33;
        h2 = h2.wrapping_mul(0xff51_afd7_ed55_8ccd);
        h2 ^= h2 >> 33;
        [h1, h2]
    }

    /// Walks the decomposition table's units once and resolves the gating
    /// for the given presence bitsets and degree-support masks: the flat
    /// per-entry ticket list (orientation bits packed) plus, in the same
    /// walk, the per-word transposition of every active contribution into
    /// its target word's write class. Shared by both prologue variants.
    ///
    /// `entries` (public i/j) drives the p-run walk's degree lookups —
    /// degrees are layout-independent — while `presence_entries` must
    /// carry the i/j relabeling that matches the presence bitsets' index
    /// space (public entries for the public prologue, the class-order's
    /// relabeled table for the class one). The two tables share order,
    /// index space, and `decomp_start`s, so the ticket's entry index is
    /// layout-independent either way. `decomp_tbl`/`pos_degree` are the
    /// working layout's scatter indices and position→degree map, and
    /// `space` is the result buffer's position count.
    fn build_gating(
        table: &FeasibleDecompositions<U>,
        entries: &[Entry],
        presence_entries: &[Entry],
        decomp_tbl: &[u32],
        a_present: &[u64],
        b_present: &[u64],
        a_deg: [u64; 2],
        b_deg: [u64; 2],
    ) -> KernelGating {
        debug_assert!(
            table.len() < (1 << 30),
            "entry indices must fit a ticket's 30 bits"
        );
        let mut tickets: Vec<u32> = Vec::new();
        let mut unit_words: Vec<u32> = Vec::new();
        let mut units: Vec<UnitRange> = Vec::with_capacity(table.units().len());
        // Transient per-unit rel bitset over the unit's degree slice,
        // cleared by iterating its set bits after extraction.
        let mut rel_bits: Vec<u64> = Vec::new();
        for unit in table.units().iter() {
            let t = unit.target as usize;
            let rs = table.degree_start(t) as u32;
            let ticket_start = tickets.len() as u32;
            let mut pairs = 0u32;
            rel_bits.clear();
            rel_bits.resize((table.degree_start(t + 1) as u32 - rs).div_ceil(64) as usize, 0);
            let mut cur_p = u8::MAX;
            let mut run_start = unit.start;
            // Real entries only: `unit.end` is the trailing sentinel's
            // slot (its decomp_start closes the last run's last
            // decomposition range via the +1 span).
            for ei in unit.start..unit.end {
                let p = table.degree_of(entries[ei as usize].i as usize) as u8;
                if p == cur_p {
                    continue;
                }
                if cur_p != u8::MAX {
                    Self::push_run(
                        presence_entries,
                        decomp_tbl,
                        a_present,
                        b_present,
                        &mut tickets,
                        &mut rel_bits,
                        &mut pairs,
                        a_deg,
                        b_deg,
                        cur_p,
                        t,
                        run_start,
                        ei,
                    );
                }
                cur_p = p;
                run_start = ei;
            }
            if cur_p != u8::MAX {
                Self::push_run(
                    presence_entries,
                    decomp_tbl,
                    a_present,
                    b_present,
                    &mut tickets,
                    &mut rel_bits,
                    &mut pairs,
                    a_deg,
                    b_deg,
                    cur_p,
                    t,
                    run_start,
                    unit.end,
                );
            }
            // Extract the unit's ascending word positions (pos = rs + rel)
            // and clear the bitset for the next unit.
            let mut emitted = 0u32;
            for (w, bits) in rel_bits.iter_mut().enumerate() {
                let mut b = *bits;
                while b != 0 {
                    let bit = b.trailing_zeros();
                    b &= b - 1;
                    unit_words.push(rs + (w as u32) * 64 + bit);
                    emitted += 1;
                }
                *bits = 0;
            }
            // A unit with no active contributions (every run gated out, or
            // every active entry's row empty) is omitted: no tickets, no
            // words, no work.
            if emitted > 0 {
                units.push(UnitRange {
                    rs,
                    td: unit.target,
                    ticket_start,
                    ticket_end: tickets.len() as u32,
                    pairs,
                });
            } else {
                debug_assert_eq!(pairs, 0);
                tickets.truncate(ticket_start as usize);
            }
        }
        // Units of one degree are ordered by CONTENT BYTES (the table's
        // unit order), not by position — the per-unit concatenation is not
        // ascending. Sort the union so the flat list is globally ascending:
        // the sink/scatter-set walk order must stay byte-identical to the
        // per-word fan-in's (which emitted active positions in ascending
        // order).
        unit_words.sort_unstable();
        let total_pairs = units.iter().map(|u| u.pairs as usize).sum();
        KernelGating {
            units: Arc::from(units),
            tickets: Arc::from(tickets),
            unit_words: Arc::from(unit_words),
            total_pairs,
        }
    }

    /// Class-space variant of [`Self::kernel_prologue_cached`]: the
    /// support lists are class-indexed, so the presence bitsets are
    /// class-positioned, the degree masks read through the ordering's
    /// relabeled degree table, and the per-word transposition targets
    /// class positions. The memo key — the exact class-indexed support
    /// lists — is the class-space image of the public variant's, so each
    /// working mode resolves its own gating entries (the class order is a
    /// pure function of the basis, so the two modes' transpositions are
    /// the same permutation apart).
    fn kernel_prologue_cached_class(
        a_series: &LieSeries<T, U>,
        a_nonzero_cls: &[usize],
        b_nonzero_cls: &[usize],
        order: &ClassOrder,
        cache: &mut GatingCache,
    ) -> KernelGating {
        Self::kernel_prologue_cached_class_core(
            &a_series.feasible_decompositions,
            a_series.basis.len(),
            a_nonzero_cls,
            b_nonzero_cls,
            order,
            cache,
        )
    }

    /// T-free core of [`Self::kernel_prologue_cached_class`] (see
    /// `sweep_words_serial_core` for why the letter type must not reach
    /// the parallel kernels).
    fn kernel_prologue_cached_class_core(
        table: &FeasibleDecompositions<U>,
        basis_len: usize,
        a_nonzero_cls: &[usize],
        b_nonzero_cls: &[usize],
        order: &ClassOrder,
        cache: &mut GatingCache,
    ) -> KernelGating {

        // Full-support fast path (class space): same predicate and cutoff
        // logic as the public prologue — the supports must be EXACTLY the
        // prefix `0..cutoff` (length alone misfires on batch-recorded node
        // lists containing degree-`max_degree` positions; see the public
        // path's comment and the pooled-dag leaf-reuse regression test),
        // per-unit gating cached on the ordering.
        let cutoff = table.degree_start(table.max_degree());
        let full_prefix = |s: &[usize]| s.len() == cutoff && s.iter().all(|&p| p < cutoff);
        if full_prefix(a_nonzero_cls) && full_prefix(b_nonzero_cls) {
            let (units, tickets, unit_words) = order.full_support_gating_class(table);
            return KernelGating {
                total_pairs: units.iter().map(|u| u.pairs as usize).sum(),
                units,
                tickets,
                unit_words,
            };
        }

        let key = (
            Self::support_fingerprint(a_nonzero_cls),
            Self::support_fingerprint(b_nonzero_cls),
        );
        if let Some((units, tickets, unit_words, total)) = cache.get(key) {
            return KernelGating {
                units: units.clone(),
                tickets: tickets.clone(),
                unit_words: unit_words.clone(),
                total_pairs: *total,
            };
        }
        // Class-positioned presence bitsets (indexed by class positions)
        // and relabeled degrees; the rest of the gating walk is shared
        // with the public prologue (see `build_gating`).
        let words = basis_len.div_ceil(64);
        let mut presence = vec![0u64; 2 * words];
        let (a_present, b_present) = presence.split_at_mut(words);
        let mut a_deg = [0u64; 2];
        let mut b_deg = [0u64; 2];
        for &i in a_nonzero_cls {
            a_present[i / 64] |= 1u64 << (i % 64);
            let d = order.degree_cls()[i] as usize;
            debug_assert!(d < 128, "degree masks cover degrees 0..127");
            a_deg[d / 64] |= 1u64 << (d % 64);
        }
        for &j in b_nonzero_cls {
            b_present[j / 64] |= 1u64 << (j % 64);
            let d = order.degree_cls()[j] as usize;
            debug_assert!(d < 128, "degree masks cover degrees 0..127");
            b_deg[d / 64] |= 1u64 << (d % 64);
        }
        let value = Self::build_gating(
            table,
            table.entries(),
            order.entries_cls(),
            order.decomp_cls(),
            a_present,
            b_present,
            a_deg,
            b_deg,
        );
        let KernelGating {
            ref units,
            ref tickets,
            ref unit_words,
            total_pairs,
        } = value;
        cache.insert(
            key,
            (
                Arc::clone(units),
                Arc::clone(tickets),
                Arc::clone(unit_words),
                total_pairs,
            ),
        );
        value
    }
    /// Resolves one maximal p-run of a unit: run-level gating on the
    /// degree-support masks, then the per-entry presence tests whose
    /// results become the run's packed tickets — and, in the same walk,
    /// each active entry's decomposition row marks its target rels in the
    /// owning unit's transient bitset and counts the unit's pairs.
    /// Tickets are built in table order, so each word's contribution
    /// sequence is a subsequence of the entry stream in table order:
    /// per-word float summation order is provably unchanged.
    #[allow(clippy::too_many_arguments)]
    fn push_run(
        presence_entries: &[Entry],
        decomp_tbl: &[u32],
        a_present: &[u64],
        b_present: &[u64],
        tickets: &mut Vec<u32>,
        rel_bits: &mut [u64],
        pairs: &mut u32,
        a_deg: [u64; 2],
        b_deg: [u64; 2],
        p: u8,
        t: usize,
        run_start: u32,
        run_end: u32,
    ) {
        let pu = p as usize;
        let qu = t - pu;
        let o1 = a_deg[pu / 64] >> (pu % 64) & 1 != 0 && b_deg[qu / 64] >> (qu % 64) & 1 != 0;
        let o2 = b_deg[pu / 64] >> (pu % 64) & 1 != 0 && a_deg[qu / 64] >> (qu % 64) & 1 != 0;
        if !o1 && !o2 {
            return;
        }
        // The same bitmap expressions the sweeps used to re-evaluate per
        // kernel call, resolved once here: `o1`/`o2` ANDed in identically,
        // presence tested per index.
        for ei in run_start..run_end {
            // i/j in the presence bitsets' index space (public for the
            // public prologue, class positions for the class one).
            let entry = &presence_entries[ei as usize];
            let (i, j) = (entry.i as usize, entry.j as usize);
            let p_active = o1
                && a_present[i / 64] & (1u64 << (i % 64)) != 0
                && b_present[j / 64] & (1u64 << (j % 64)) != 0;
            let q_active = o2
                && a_present[j / 64] & (1u64 << (j % 64)) != 0
                && b_present[i / 64] & (1u64 << (i % 64)) != 0;
            if !p_active && !q_active {
                continue;
            }
            tickets.push(
                ei | if p_active { TICKET_P_ACTIVE } else { 0 }
                    | if q_active { TICKET_Q_ACTIVE } else { 0 },
            );
            // Mark the entry's row targets in the owning unit's bitset
            // (rels are degree-slice relative — the same in both working
            // spaces) and count the pairs.
            let from = entry.decomp_start as usize;
            let to = presence_entries[ei as usize + 1].decomp_start as usize;
            for &rel in &decomp_tbl[from..to] {
                rel_bits[(rel / 64) as usize] |= 1u64 << (rel % 64);
                *pairs += 1;
            }
        }
    }
    fn sweep_words_serial<S: ScatterSink>(
        a_series: &LieSeries<T, U>,
        a_coefficients: &[U],
        b_coefficients: &[U],
        gating: &KernelGating,
        result_coefficients: &mut [U],
        entries: &[Entry],
        decomp_rels: &[u32],
        sink: &mut S,
    ) {
        Self::sweep_words_serial_core(
            a_series.feasible_decompositions.decomp_coeffs(),
            a_coefficients,
            b_coefficients,
            gating,
            result_coefficients,
            entries,
            decomp_rels,
            sink,
        );
    }

    /// T-free core of [`Self::sweep_words_serial`]: everything the sweep
    /// touches is basis-derived integer tables plus U values — the letter
    /// type never enters the kernel (this is what lets the collecting
    /// rebuild run level-parallel; see `CommutatorDag::ensure_lists_steady`).
    fn sweep_words_serial_core<S: ScatterSink>(
        decomp_coefficients: &[U],
        a_coefficients: &[U],
        b_coefficients: &[U],
        gating: &KernelGating,
        result_coefficients: &mut [U],
        entries: &[Entry],
        decomp_rels: &[u32],
        sink: &mut S,
    ) {
        // Single-phase per-unit sweep: each unit's active entries are
        // visited in table order, each term computed exactly once and its
        // row contributions added straight into the result buffer. A
        // word's adds happen inside its one owning unit, in table-entry
        // order — exactly the serial scatter's per-word `+=` sequence —
        // and callers' buffers are accumulated into (their current value
        // starts the sum). Degree-slice starts are identical in both
        // working layouts, so `rs + rel` addresses the working-space
        // result directly.
        for u in gating.units.iter() {
            for tp in u.ticket_start as usize..u.ticket_end as usize {
                let ticket = gating.tickets[tp];
                let e = (ticket & TICKET_INDEX_MASK) as usize;
                let entry = entries[e];
                let p_active = ticket & TICKET_P_ACTIVE != 0;
                let q_active = ticket & TICKET_Q_ACTIVE != 0;
                let (i, j) = (entry.i as usize, entry.j as usize);
                // Orientation (`a = j, b = i` only) negates: `[basis[j],
                // basis[i]]` is the negation of the stored decomposition.
                // `raw_mul`/`raw_add_assign` skip the float wrappers'
                // per-op NaN checks (raw-float fast path, see `raw_mul`'s
                // NaN policy).
                let term = if p_active {
                    let mut t = raw_mul(&a_coefficients[i], &b_coefficients[j]);
                    if q_active {
                        raw_add_assign(&mut t, &-raw_mul(&a_coefficients[j], &b_coefficients[i]));
                    }
                    t
                } else {
                    -raw_mul(&a_coefficients[j], &b_coefficients[i])
                };
                let from = entry.decomp_start as usize;
                let to = entries[e + 1].decomp_start as usize;
                let base = u.rs as usize;
                for (k, &rel) in decomp_rels[from..to].iter().enumerate() {
                    let pos = base + rel as usize;
                    raw_add_assign(
                        &mut result_coefficients[pos],
                        &raw_mul(&decomp_coefficients[from + k], &term),
                    );
                }
            }
        }
        // All units are complete: report every active position (globally
        // ascending — the same sequence the per-word fan-in reported). No
        // unit ever writes another's words.
        for &pos in gating.unit_words.iter() {
            sink.scatter(pos as usize);
        }
    }

    pub fn commutator_assign(&mut self, other: &Self)
    where
        U: Send + Sync,
    {
        let mut coefficients = vec![U::default(); self.coefficients.len()];
        LieSeries::commutator_coefficients(
            self,
            &self.coefficients,
            &other.coefficients,
            &mut coefficients,
        );
        self.coefficients = coefficients;
    }
}

impl<
    T: Clone + Ord + Generator + Hash,
    U: Clone
        + Default
        + One
        + Zero
        + Eq
        + MulAssign
        + Neg<Output = U>
        + Hash
        + AddAssign
        + Send
        + Sync
        + 'static,
> Commutator<&Self> for LieSeries<T, U>
{
    type Output = Self;

    /// Calculates the lie bracket `[A, B]` for a lie series for terms within the commutator basis.
    fn commutator(&self, other: &Self) -> Self::Output {
        let mut coefficients = vec![U::default(); self.coefficients.len()];
        LieSeries::<T, U>::commutator_coefficients(
            self,
            // other,
            &self.coefficients,
            &other.coefficients,
            &mut coefficients,
        );
        self.with_coefficients(coefficients)
    }
}

#[cfg(test)]
mod test {
    use ordered_float::NotNan;
    use rstest::rstest;

    /// The work-adaptive slot policy must reproduce the measured regimes:
    /// tiny folds walk serially at any pool size, mid folds land on their
    /// measured sweet spot, wide folds fill the pool, and the pool/pack
    /// caps still bind.
    #[test]
    fn work_adaptive_slots_matches_measured_regimes() {
        use super::work_adaptive_slots;
        // 12x2 per-fold (66 entries/fold): serial at every pool size.
        assert_eq!(work_adaptive_slots(32, 32, 66), 1);
        assert_eq!(work_adaptive_slots(2, 32, 66), 1);
        // 8x3 (~700 entries/fold, + block elements in the batch path):
        // 1-2 slots.
        assert_eq!(work_adaptive_slots(32, 32, 700), 1);
        assert!(work_adaptive_slots(32, 32, 2000) <= 2, "8x3 regime: 2000 entries must fund at most 2 slots");
        // 3x8 (~74K planned entries/fold): the QUANTUM lands it near the
        // measured 8-16t sweet spot — 19 slots at a 32t pool (74_000/3750).
        assert_eq!(work_adaptive_slots(32, 32, 74_000), 19);
        // ... and stays in that neighborhood across the regime's ticket
        // variance (the real 3x8 fold's planned count wobbles around 74K).
        let s = work_adaptive_slots(32, 32, 70_000);
        assert!((16..=21).contains(&s), "3x8 regime: 70K entries must land in 16..=21 slots, got {s}");
        // 2x12 (~257K entries/fold): the full pool.
        assert_eq!(work_adaptive_slots(32, 32, 256_858), 32);
        // Pool and pack caps bind before the work term.
        assert_eq!(work_adaptive_slots(8, 3, 256_858), 3);
        assert_eq!(work_adaptive_slots(1, 32, 256_858), 1);
        assert_eq!(work_adaptive_slots(32, 1, 256_858), 1);
        // Degenerate: no planned work → serial.
        assert_eq!(work_adaptive_slots(32, 32, 0), 1);
    }

    /// INVARIANT CHECK on the real tables: within one degree-target slice,
    /// two different units must never scatter onto the same basis word (the
    /// batch's packs cut at unit boundaries and rely on per-word
    /// single-writer accumulation).
    #[test]
    fn debug_real_table_unit_word_ownership() {
        for (d, m) in [
            (2usize, 4usize),
            (3usize, 4usize),
            (3usize, 5usize),
            (3usize, 8usize),
        ] {
            let basis = lyndon_rs::lyndon::LyndonBasis::<u8>::new(
                d,
                lyndon_rs::lyndon::Sort::Lexicographical,
            )
            .generate_basis(m);
            let len = basis.len();
            let a = LieSeries::<u8, NotNan<f64>>::new(
                basis.clone(),
                (0..len)
                    .map(|k| NotNan::new((k % 17) as f64 / 19.0 - 0.4).unwrap())
                    .collect(),
            );
            let _b = LieSeries::<u8, NotNan<f64>>::new(
                basis,
                (0..len)
                    .map(|k| NotNan::new((k % 13) as f64 / 17.0 - 0.3).unwrap())
                    .collect(),
            );
            let fd = &a.feasible_decompositions;
            let mut owner: std::collections::HashMap<(u8, u32), usize> =
                std::collections::HashMap::new();
            let mut violations = 0usize;
            for (uid, unit) in fd.units().iter().enumerate() {
                let span = fd.entry_span(unit);
                for (entry, next) in span[..span.len() - 1].iter().zip(span[1..].iter()) {
                    let from = entry.decomp_start as usize;
                    let to = next.decomp_start as usize;
                    for &rel in &fd.decomp_indices_rel()[from..to] {
                        let key = (unit.target, rel);
                        match owner.get(&key) {
                            Some(prev) if *prev != uid => {
                                violations += 1;
                                if violations <= 6 {
                                    println!(
                                        "VIOLATION d={d} m={m}: degree {} rel {} by units {} and {}",
                                        unit.target, rel, prev, uid
                                    );
                                }
                            }
                            _ => {
                                owner.insert(key, uid);
                            }
                        }
                    }
                }
            }
            println!(
                "d={d} m={m}: units={} rels={} violations={}",
                fd.units().len(),
                fd.decomp_indices_rel().len(),
                violations
            );
            assert_eq!(violations, 0);
        }
    }

    /// Anagram-class scatter locality of the commutation kernel's write set:
    /// per anagram unit, how many distinct 64-byte lines its stores touch and
    /// how far they spread inside the public degree slice, versus the same
    /// stores routed into the class-contiguous (internal) layout.
    #[test]
    #[ignore = "layout stats: run explicitly with --ignored --nocapture"]
    fn dump_scatter_locality_stats() {
        for (d, m) in [(2usize, 12usize), (3, 8), (4, 10), (12, 2)] {
            let basis = lyndon_rs::LyndonBasis::<u8>::new(d, lyndon_rs::Sort::Lexicographical);
            let words = basis.generate_basis(m);
            let len = words.len();
            let a = LieSeries::new(words.clone(), vec![NotNan::new(1.0).unwrap(); len]);
            let table = &a.feasible_decompositions;
            let decomp = table.decomp_indices_rel();
            let entries = table.entries();
            let mut cur_lines = 0usize;
            let mut cls_lines = 0usize;
            let mut total_rmw = 0usize;
            let mut spread_max = 0usize;
            let mut report = String::new();
            for unit in table.units() {
                let (rs, re) = table.degree_range(unit.target);
                let slice = (re - rs) as usize;
                let mut words_seen = std::collections::BTreeSet::new();
                let mut rmw = 0usize;
                for ei in unit.start..unit.end {
                    let from = entries[ei as usize].decomp_start as usize;
                    let to = entries[ei as usize + 1].decomp_start as usize;
                    for &rel in &decomp[from..to] {
                        words_seen.insert(rel as usize);
                    }
                    rmw += to - from;
                }
                let distinct = words_seen.len();
                if distinct == 0 {
                    continue;
                }
                let lo = *words_seen.first().unwrap();
                let hi = *words_seen.last().unwrap();
                let spread = hi - lo + 1;
                let lines_touched = words_seen
                    .iter()
                    .map(|&w| (rs as usize + w) / 8)
                    .collect::<std::collections::HashSet<_>>()
                    .len();
                let lines_cls = (distinct + 7) / 8;
                cur_lines += lines_touched;
                cls_lines += lines_cls;
                total_rmw += rmw;
                spread_max = spread_max.max(spread);
                report.push_str(&format!(
                    "    unit t={:2} words={:5} rmw={:6} slice={:7} spread={:6} ({:5.1}% of slice) lines_cur={:5} lines_cls={:5}\n",
                    unit.target,
                    distinct,
                    rmw,
                    slice,
                    spread,
                    100.0 * spread as f64 / slice.max(1) as f64,
                    lines_touched,
                    lines_cls,
                ));
            }
            println!(
                "{d}x{m}: len={len} store_rmw={total_rmw} lines_cur_sum={cur_lines} lines_cls_sum={cls_lines} ({:.2}x reduction) max_unit_spread={spread_max}",
                cur_lines as f64 / cls_lines.max(1) as f64
            );
            print!("{report}");
        }
    }

    use super::*;

    #[rstest]
    #[case(2, 2, vec![1, 2, 3], vec![4, 5, 6], vec![0, 0, -3])]
    #[case(2, 2, vec![3, 2, 1], vec![1, 2, 3], vec![0, 0, 4])]
    #[case(2, 3, vec![1, 2, 3, 4, 5], vec![6, 7, 8, 9, 10], vec![0, 0, -5, -10, 5])]
    #[case(3, 3,
        vec![1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        vec![5, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        vec![0, 0, 0, -7, -14, -7, 0, 0, 0, 0, 0, 0, 0, 0])]
    #[case(3, 3,
        vec![1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        vec![0, 0, 0, -7, -14, -7, 0, 0, 0, 0, 0, 0, 0, 0],
        vec![0, 0, 0, 0, 0, 0, -7, -14, 14, 14, 49, 42, -14, 21])]
    #[case(3, 4,
        vec![
            1, 2, 3, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0
        ],
        vec![
            5, 3, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0
        ],
        vec![
            0, 0, 0, -7, -14, -7, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0
        ],
    )]
    #[case(3, 4,
        vec![
            1, 2, 3, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0
        ],
        vec![
            0, 0, 0, -7, -14, -7, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0
        ],
        vec![
            0, 0, 0, 0, 0, 0, -7, -14,
            14, 14, 49, 42, -14, 21, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
        ],
    )]
    fn test_lie_series_commutation(
        #[case] num_generators: usize,
        #[case] basis_depth: usize,
        #[case] a_coefficients: Vec<i128>,
        #[case] b_coefficients: Vec<i128>,
        #[case] expected_coefficients: Vec<i128>,
    ) {
        use lyndon_rs::lyndon::{LyndonBasis, Sort};
        use num_rational::Ratio;

        let basis = LyndonBasis::<u8>::new(num_generators, Sort::Lexicographical)
            .generate_basis(basis_depth);
        let a_coefficients = a_coefficients
            .into_iter()
            .map(Ratio::<i128>::from_integer)
            .collect::<Vec<_>>();
        let a = LieSeries::new(basis.clone(), a_coefficients);

        let b_coefficients = b_coefficients
            .into_iter()
            .map(Ratio::<i128>::from_integer)
            .collect::<Vec<_>>();
        let b = LieSeries::new(basis.clone(), b_coefficients);
        let expected_coefficients = expected_coefficients
            .into_iter()
            .map(Ratio::<i128>::from_integer)
            .collect::<Vec<_>>();

        let series = comm![a, b];
        assert_eq!(series.coefficients.len(), expected_coefficients.len());
        dbg!(&series.coefficients);
        for (i, c) in series.coefficients.iter().enumerate() {
            assert_eq!(
                *c, expected_coefficients[i],
                "{i}: {c:?} != {:?}",
                expected_coefficients[i]
            );
        }
    }

    #[test]
    fn decompositions_contain_no_zero_coefficients() {
        use lyndon_rs::lyndon::{LyndonBasis, Sort};
        use ordered_float::NotNan;

        for (d, m) in [(2usize, 6usize), (3, 5), (4, 4)] {
            let words = LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
            let series = LieSeries::<u8, NotNan<f64>>::new(words, Vec::new());
            assert!(
                series
                    .feasible_decompositions
                    .iter_entries()
                    .flat_map(|(_, _, _, coefficients)| coefficients)
                    .all(|c| !c.is_zero()),
                "zero coefficient found in decompositions for d={d}, m={m}"
            );
        }
    }

    /// Independently recomputes the decomposition of every feasible canonical
    /// pair and pins the compact table's contents (indices, coefficients,
    /// ordering, feasibility, and canonicality) against it, plus the mirrored
    /// (negated) query path.
    #[test]
    fn feasible_decompositions_match_independent_recomputation() {
        use lyndon_rs::lyndon::{LyndonBasis, Sort};
        use ordered_float::NotNan;

        for (d, m) in [(2usize, 6usize), (3, 5), (4, 4)] {
            let words = LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
            let series = LieSeries::<u8, NotNan<f64>>::new(words, Vec::new());
            let degree = |x: usize| series.commutator_basis[x].degree();
            let index_of = |term: &CommutatorTerm<NotNan<f64>, u8>| {
                series.commutator_basis_index_map[&term.unit_hash()]
            };
            let basis_set: HashSet<_> = series
                .commutator_basis
                .iter()
                .map(CommutatorTerm::unit_hash)
                .collect();

            let mut feasible = 0;
            for unit in series.feasible_decompositions.units() {
                let (p, _q, target) = (unit.p_mask, unit.target as usize, unit.target as usize);
                assert_eq!(
                    p.iter().map(|w| w.count_ones()).sum::<u32>() as usize >= 1,
                    true,
                    "empty unit"
                );
                assert_eq!(
                    unit.gamma.iter().map(|&x| x as usize).sum::<usize>(),
                    target,
                    "unit target degree mismatch"
                );
                let entries: Vec<_> = series
                    .feasible_decompositions
                    .unit_iter(unit)
                    .map(|(e, _, _)| (e.i as usize, e.j as usize))
                    .collect();
                assert!(
                    entries.windows(2).all(|w| w[0] < w[1]),
                    "unit {:?} entries not sorted by (i, j)",
                    unit.gamma
                );

                for (i, j) in entries {
                    assert!(i < j, "non-canonical pair stored");
                    assert!(
                        p[(degree(i) / 64) as usize] >> (degree(i) % 64) & 1 == 1,
                        "pair ({i}, {j}) in a unit of the wrong degrees"
                    );
                    assert_eq!(
                        degree(i) + degree(j),
                        target,
                        "pair ({i}, {j}) in a unit of the wrong target degree"
                    );
                    assert!(degree(i) + degree(j) <= m, "infeasible pair stored");

                    let mut term = comm![&series.commutator_basis[i], &series.commutator_basis[j]];
                    term.lyndon_sort();
                    let expected: Vec<_> = term
                        .lyndon_basis_decomposition(&basis_set)
                        .into_iter()
                        .filter(|x| !x.coefficient().is_zero())
                        .collect();

                    let expected_indices: Vec<u32> =
                        expected.iter().map(|x| index_of(x) as u32).collect();
                    let expected_coefficients: Vec<_> =
                        expected.iter().map(|x| x.coefficient().clone()).collect();

                    // The canonical query returns the stored data unflagged;
                    // comparing through `decomposition` exercises the
                    // degree-block lookup as well as the contents.
                    let (canonical_indices, canonical_coefficients, swapped) =
                        series.decomposition(i, j).expect("canonical pair query");
                    assert!(!swapped, "(i={i}, j={j}) canonical query flagged");
                    assert_eq!(
                        canonical_indices,
                        &expected_indices[..],
                        "(i={i}, j={j}) indices"
                    );
                    assert_eq!(
                        canonical_coefficients,
                        &expected_coefficients[..],
                        "(i={i}, j={j}) coeffs"
                    );
                    feasible += 1;

                    // The mirrored query returns the same (canonical) data,
                    // flagged so the caller negates it into
                    // [basis[j], basis[i]] orientation.
                    let (mirrored_indices, mirrored_coefficients, swapped) =
                        series.decomposition(j, i).expect("mirrored pair query");
                    assert!(swapped, "(i={i}, j={j}) mirrored query not flagged");
                    assert_eq!(mirrored_indices, &expected_indices[..]);
                    assert_eq!(
                        mirrored_coefficients,
                        &expected_coefficients[..],
                        "(i={i}, j={j}) mirrored data differs from canonical"
                    );
                }
            }
            assert_eq!(series.feasible_decompositions.len(), feasible);
        }
    }

    /// The fused kernel evaluates both bracket orientations of every canonical
    /// pair, so `[A, A]` must vanish exactly — every pair contributes
    /// `c * (a_min * a_max - a_max * a_min) = 0`.
    #[test]
    fn commutator_of_series_with_itself_is_zero() {
        use lyndon_rs::lyndon::{LyndonBasis, Sort};
        use ordered_float::NotNan;

        for (d, m) in [(2usize, 6usize), (3, 5)] {
            let words: Vec<LyndonWord<u8>> =
                LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
            let coefficients: Vec<_> = (0..words.len())
                .map(|i: usize| NotNan::new(((i % 7 + 1) as f64) / 3.0 - 1.1).unwrap())
                .collect();
            let series: LieSeries<u8, NotNan<f64>> = LieSeries::new(words, coefficients);
            let result: LieSeries<u8, NotNan<f64>> = series.commutator(&series);
            assert!(
                result.coefficients.iter().all(|c| c.is_zero()),
                "non-zero coefficient in [A, A] for d={d}, m={m}"
            );
        }
    }

    /// ADVERSARIAL (per-unit division): the gating's per-unit structure
    /// must be LOSSLESS and ORDERED against an independent walk of the
    /// table with presence resolved straight from the support lists:
    /// (a) the unit ticket ranges partition the flat ticket list without
    /// overlap, (b) each unit's orientation flags match the independent
    /// presence resolution per entry, (c) the unit word sets PARTITION the
    /// active word positions — pairwise disjoint, ascending per unit, and
    /// their concatenation ascending (so CollectSink order is unchanged) —
    /// and equal the set of positions the independent walk produces, (d)
    /// each unit's contribution sequence — its tickets expanded in order,
    /// each entry's row (rel, dp) pairs — is exactly the independent
    /// walk's subsequence for that unit's words, in table-entry order (the
    /// bit-exactness contract), and (e) total_pairs matches. A
    /// misbucketed, lost, duplicated, or reordered contribution fails
    /// here with a distinct message per invariant.
    #[test]
    fn write_class_gating_transposition_is_lossless_and_ordered() {
        use lyndon_rs::lyndon::{LyndonBasis, Sort};
        use ordered_float::NotNan;

        use super::ClassOrderedCommutation;

        for (d, m) in [(2usize, 4usize), (3usize, 5usize), (3usize, 8usize)] {
            let words: Vec<LyndonWord<u8>> =
                LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
            let len = words.len();
            let series: LieSeries<u8, NotNan<f64>> = LieSeries::new(
                words.clone(),
                (0..len)
                    .map(|i| NotNan::new(((i % 7 + 1) as f64) / 5.0 - 0.7).unwrap())
                    .collect(),
            );
            let table = &series.feasible_decompositions;
            let entries = table.entries();
            let decomp = table.decomp_indices_rel();
            let degrees = table.index_degrees();
            let order = series.class_order();

            // Independent support pairs: full, one-sided, random halves,
            // single letters, empty.
            let cutoff = table.degree_start(table.max_degree());
            let mut support_sets: Vec<Vec<usize>> = vec![Vec::new()];
            support_sets.push((0..cutoff).collect());
            support_sets.push((0..cutoff).filter(|i| i % 2 == 0).collect());
            support_sets.push((0..cutoff).filter(|i| i % 3 != 1).collect());
            support_sets.push((0..cutoff).filter(|i| degrees[*i] == 1).collect());
            support_sets.push(vec![0]);

            for a_support in &support_sets {
                for b_support in &support_sets {
                    // Independent active-contribution walk over the table,
                    // in table-entry order: presence straight from the
                    // support lists (degree-mask gating is implied by
                    // membership).
                    let present = |s: &[usize]| {
                        let mut p = vec![false; len];
                        for &i in s {
                            p[i] = true;
                        }
                        p
                    };
                    let (ap, bp) = (present(a_support), present(b_support));
                    // (pos, entry, dp) in table-entry order.
                    let mut expected: Vec<(usize, usize, usize)> = Vec::new();
                    for (ei, entry) in entries[..entries.len() - 1].iter().enumerate() {
                        let (i, j) = (entry.i as usize, entry.j as usize);
                        let p_active = ap[i] && bp[j];
                        let q_active = ap[j] && bp[i];
                        if !p_active && !q_active {
                            continue;
                        }
                        let from = entry.decomp_start as usize;
                        let to = entries[ei + 1].decomp_start as usize;
                        let rs = table.degree_start(degrees[i] as usize + degrees[j] as usize);
                        for (k, &rel) in decomp[from..to].iter().enumerate() {
                            expected.push((rs + rel as usize, ei, from + k));
                        }
                    }
                    // Active entries with orientation, table order (the
                    // tickets must be exactly this list).
                    let mut expected_tickets: Vec<(usize, bool, bool)> = Vec::new();
                    for (ei, entry) in entries[..entries.len() - 1].iter().enumerate() {
                        let (i, j) = (entry.i as usize, entry.j as usize);
                        let p_active = ap[i] && bp[j];
                        let q_active = ap[j] && bp[i];
                        if p_active || q_active {
                            expected_tickets.push((ei, p_active, q_active));
                        }
                    }

                    // Public-mode gating: positions are public basis
                    // indices.
                    let mut cache = GatingCache::default();
                    let gating = LieSeries::<u8, NotNan<f64>>::kernel_prologue_cached(
                        &series,
                        a_support,
                        b_support,
                        &mut cache,
                    );
                    let ctx = "public";

                    // (a) Ticket-range partition: consecutive, ordered,
                    // non-overlapping, covering the ticket list.
                    let mut cursor = 0usize;
                    for (ui, u) in gating.units.iter().enumerate() {
                        assert_eq!(
                            u.ticket_start as usize, cursor,
                            "{ctx}: unit {ui} ticket range not contiguous"
                        );
                        assert!(
                            u.ticket_end > u.ticket_start,
                            "{ctx}: unit {ui} emitted with no tickets"
                        );
                        assert_eq!(
                            degrees[u.rs as usize], u.td,
                            "{ctx}: unit {ui} rs/td degree mismatch"
                        );
                        cursor = u.ticket_end as usize;
                    }
                    assert_eq!(
                        cursor,
                        gating.tickets.len(),
                        "{ctx}: ticket ranges do not cover the ticket list"
                    );

                    // (b) Orientation flags match the independent walk,
                    // ticket list = expected active entries in table order.
                    assert_eq!(
                        gating.tickets.len(),
                        expected_tickets.len(),
                        "{ctx}: active-entry count mismatch"
                    );
                    for (tp, &(ei, want_p, want_q)) in
                        expected_tickets.iter().enumerate()
                    {
                        let ticket = gating.tickets[tp];
                        assert_eq!(
                            (ticket & TICKET_INDEX_MASK) as usize, ei,
                            "{ctx}: ticket {tp} entry mismatch"
                        );
                        assert_eq!(
                            ticket & TICKET_P_ACTIVE != 0, want_p,
                            "{ctx}: p_active flag mismatch at entry {ei}"
                        );
                        assert_eq!(
                            ticket & TICKET_Q_ACTIVE != 0, want_q,
                            "{ctx}: q_active flag mismatch at entry {ei}"
                        );
                    }

                    // (c)+(d) Word-set partition AND per-unit sequence:
                    // reconstruct each unit's word set and contribution
                    // sequence from its rows; sets must be pairwise
                    // disjoint, within the unit's degree slice, and their
                    // union must equal both the gating's flat `unit_words`
                    // list and the independent walk's positions; each
                    // unit's (entry, dp) sequence must be exactly the
                    // independent walk's subsequence for the unit's words,
                    // in table-entry order (the bit-exactness contract).
                    let mut all_positions: Vec<usize> = Vec::new();
                    let mut seen_pairs: Vec<(usize, usize)> = Vec::new();
                    for (ui, u) in gating.units.iter().enumerate() {
                        let mut set: Vec<usize> = Vec::new();
                        let mut unit_pairs: Vec<(usize, usize)> = Vec::new();
                        for tp in u.ticket_start as usize..u.ticket_end as usize {
                            let ei = (gating.tickets[tp] & TICKET_INDEX_MASK) as usize;
                            let entry = entries[ei];
                            let from = entry.decomp_start as usize;
                            let to = entries[ei + 1].decomp_start as usize;
                            for (k, &rel) in decomp[from..to].iter().enumerate() {
                                let pos = u.rs as usize + rel as usize;
                                set.push(pos);
                                unit_pairs.push((ei, from + k));
                            }
                        }
                        set.sort_unstable();
                        set.dedup();
                        assert!(
                            set.iter().all(|&p| {
                                degrees[p] == u.td
                                    && (u.rs as usize..table.degree_start(u.td as usize + 1))
                                        .contains(&p)
                            }),
                            "{ctx}: unit {ui} word set outside its degree slice"
                        );
                        for &p in &set {
                            assert!(
                                !all_positions.contains(&p),
                                "{ctx}: position {p} written by two units"
                            );
                        }
                        all_positions.extend_from_slice(&set);
                        seen_pairs.extend_from_slice(&unit_pairs);
                        // The unit's sequence must contain only its own
                        // words, in table-entry order.
                        let mut last_ei = None::<usize>;
                        for &(ei, dp) in &unit_pairs {
                            let rel = decomp[dp];
                            assert!(
                                set.binary_search(&(u.rs as usize + rel as usize)).is_ok(),
                                "{ctx}: unit {ui} contributes to a word outside its set"
                            );
                            assert!(
                                last_ei.is_none_or(|prev| prev <= ei),
                                "{ctx}: unit {ui} contributions out of table order"
                            );
                            last_ei = Some(ei);
                        }
                    }
                    let mut sorted_all = all_positions.clone();
                    sorted_all.sort_unstable();
                    let got_flat: Vec<usize> =
                        gating.unit_words.iter().map(|&p| p as usize).collect();
                    assert_eq!(
                        sorted_all, got_flat,
                        "{ctx}: flat unit_words list differs from the union of the rows' positions"
                    );
                    let mut want_positions: Vec<usize> =
                        expected.iter().map(|&(p, _, _)| p).collect();
                    want_positions.sort_unstable();
                    want_positions.dedup();
                    assert_eq!(
                        got_flat, want_positions,
                        "{ctx}: unit word-set partition differs from the independent walk's positions"
                    );
                    // (d) global sequence: every unit's contributions
                    // concatenated = the independent walk's (entry, dp)
                    // pairs in table-entry order.
                    let mut want_pairs: Vec<(usize, usize)> = expected
                        .iter()
                        .map(|&(_, ei, dp)| (ei, dp))
                        .collect();
                    assert_eq!(
                        seen_pairs, want_pairs,
                        "{ctx}: unit contribution sequences differ from the independent walk"
                    );
                    let _ = &mut want_pairs;

                    // (e) Pair count.
                    assert_eq!(
                        gating.total_pairs,
                        expected.len(),
                        "{ctx}: contribution count mismatch (supports {a_support:?}/{b_support:?})"
                    );

                    // Class-mode gating: positions are class positions; the
                    // same invariants must hold against the class-space
                    // image of the expected list.
                    let inv = order.inv();
                    let a_cls: Vec<usize> = a_support.iter().map(|&i| inv[i] as usize).collect();
                    let b_cls: Vec<usize> = b_support.iter().map(|&i| inv[i] as usize).collect();
                    let gating_cls = LieSeries::<u8, NotNan<f64>>::kernel_prologue_cached_class(
                        &series,
                        &a_cls,
                        &b_cls,
                        &order,
                        &mut cache,
                    );
                    let degree_cls = order.degree_cls();
                    let decomp_cls = order.decomp_cls();
                    let entries_cls = order.entries_cls();
                    // Class-space image of the expected walk: map public
                    // positions through inv, keep table-entry order.
                    let mut expected_cls: Vec<(usize, usize, usize)> = Vec::new();
                    let perm = order.perm();
                    for (ei, entry) in entries_cls[..entries_cls.len() - 1].iter().enumerate() {
                        let (i, j) = (entry.i as usize, entry.j as usize);
                        // entries_cls endpoints are CLASS positions; the
                        // presence vectors are indexed by PUBLIC position
                        // (perm: class -> public).
                        let p_active = ap[perm[i] as usize] && bp[perm[j] as usize];
                        let q_active = ap[perm[j] as usize] && bp[perm[i] as usize];
                        if !p_active && !q_active {
                            continue;
                        }
                        let from = entry.decomp_start as usize;
                        let to = entries_cls[ei + 1].decomp_start as usize;
                        let rs = table.degree_start(
                            degree_cls[i] as usize + degree_cls[j] as usize,
                        );
                        for &rel in decomp_cls[from..to].iter() {
                            expected_cls.push((rs + rel as usize, ei, from as usize));
                        }
                    }
                    // (a)+(c) for class space: ticket ranges cover the
                    // list; the flat word list is globally ascending and
                    // equals the independent class-space walk's positions.
                    let mut cursor_cls = 0usize;
                    for (ui, u) in gating_cls.units.iter().enumerate() {
                        assert_eq!(
                            u.ticket_start as usize, cursor_cls,
                            "class: unit {ui} ticket range not contiguous"
                        );
                        cursor_cls = u.ticket_end as usize;
                    }
                    assert_eq!(cursor_cls, gating_cls.tickets.len(), "class: ticket ranges do not cover the list");
                    let mut want_cls: Vec<usize> =
                        expected_cls.iter().map(|&(p, _, _)| p).collect();
                    want_cls.sort_unstable();
                    want_cls.dedup();
                    let got_cls: Vec<usize> = gating_cls
                        .unit_words
                        .iter()
                        .map(|&p| p as usize)
                        .collect();
                    let mut last_cls = None::<usize>;
                    for &p in &got_cls {
                        assert!(
                            last_cls.is_none_or(|prev| prev < p),
                            "class: unit word list not globally ascending at {p}"
                        );
                        last_cls = Some(p);
                    }
                    assert_eq!(
                        got_cls, want_cls,
                        "class: unit word-set partition differs from the independent walk"
                    );
                    // (e) Pair count.
                    assert_eq!(
                        gating_cls.total_pairs,
                        expected_cls.len(),
                        "class: class-mode contribution count mismatch"
                    );
                }
            }
        }
    }

    /// ADVERSARIAL (write-access division): the kernel sweep must equal an
    /// INDEPENDENT per-word fan-in reference — built by walking the table's
    /// entries in table order with presence resolved straight from the
    /// support lists — bit for bit, across coefficient types (raw-float
    /// fast path and exact rationals), layouts (public direct and
    /// class-contiguous), thread counts, and adversarial inputs (all-zero,
    /// single-hot, planted exact cancellations). Distinct failure: this
    /// guards arithmetic and per-word accumulation order, not gating
    /// structure (the transposition test above).
    #[test]
    fn write_class_sweep_matches_independent_fanin_reference() {
        use lyndon_rs::lyndon::{LyndonBasis, Sort};
        use ordered_float::NotNan;
        use std::collections::BTreeMap;

        use super::ClassOrderedCommutation;

        // The reference: result[w] += c * term over the table's entries in
        // table order, term = (p_active ? a_i*b_j : 0) - (q_active ?
        // a_j*b_i : 0), entries with neither orientation skipped. Same add
        // sequence per word as the kernel's write classes, derived without
        // touching the gating.
        fn reference<U>(
            series: &LieSeries<u8, U>,
            a: &[U],
            a_present: &[bool],
            b: &[U],
            b_present: &[bool],
        ) -> Vec<U>
        where
            U: Clone + Default + Mul<Output = U> + AddAssign + Neg<Output = U> + PartialEq,
        {
            let table = &series.feasible_decompositions;
            let entries = table.entries();
            let coeffs = table.decomp_coeffs();
            let mut acc: BTreeMap<usize, U> = BTreeMap::new();
            for (ei, entry) in entries[..entries.len() - 1].iter().enumerate() {
                let (i, j) = (entry.i as usize, entry.j as usize);
                let p_active = a_present[i] && b_present[j];
                let q_active = a_present[j] && b_present[i];
                if !p_active && !q_active {
                    continue;
                }
                let term = if p_active {
                    let mut t = a[i].clone() * b[j].clone();
                    if q_active {
                        t += -(a[j].clone() * b[i].clone());
                    }
                    t
                } else {
                    -(a[j].clone() * b[i].clone())
                };
                let from = entry.decomp_start as usize;
                let to = entries[ei + 1].decomp_start as usize;
                let rs = table.degree_start(
                    table.degree_of(i) + table.degree_of(j),
                );
                for dp in from..to {
                    // Absolute decomp position -> the row's target: the
                    // relative array is stored per entry in row order.
                    let rel = table.decomp_indices_rel()[dp];
                    let w = rs + rel as usize;
                    let contrib = coeffs[dp].clone() * term.clone();
                    *acc.entry(w).or_insert_with(U::default) += contrib;
                }
            }
            let mut out = vec![U::default(); series.coefficients.len()];
            for (w, v) in acc {
                out[w] = v;
            }
            out
        }

        for (d, m) in [(2usize, 4usize), (3usize, 6usize)] {
            let words: Vec<LyndonWord<u8>> =
                LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
            let len = words.len();
            let table_deg = |i: usize| words[i].len();
            let _ = table_deg;

            // Adversarial value sets: zeros (support holes), signed values
            // that plant exact cancellations between entries feeding the
            // same word, and asymmetric magnitudes.
            let mk = |f: &dyn Fn(usize) -> f64| {
                (0..len)
                    .map(|i| NotNan::new(f(i)).unwrap())
                    .collect::<Vec<_>>()
            };
            let a_sets = [
                mk(&|_| 0.0),
                mk(&|i| if i == 0 { 1.0 } else { 0.0 }),
                mk(&|i| ((i % 5 + 1) as f64) / 3.0 - 1.0),
                mk(&|i| if i % 2 == 0 { 1.0 } else { -1.0 }),
            ];
            let b_sets = [
                mk(&|_| 0.0),
                mk(&|i| if i == 1 { 1.0 } else { 0.0 }),
                mk(&|i| ((i * 3 + 2) % 7) as f64 / 2.0 - 1.5),
                mk(&|i| if i % 3 == 0 { 2.0 } else { -0.5 }),
            ];

            let cutoff_len = len; // full lists; the kernel filters by value
            for a_coefficients in &a_sets {
                for b_coefficients in &b_sets {
                    let a: LieSeries<u8, NotNan<f64>> =
                        LieSeries::new(words.clone(), a_coefficients.clone());
                    let b: LieSeries<u8, NotNan<f64>> =
                        LieSeries::new(words.clone(), b_coefficients.clone());
                    let a_nonzero = a.nonzero_coefficient_indices(a_coefficients);
                    let b_nonzero = b.nonzero_coefficient_indices(b_coefficients);
                    let ap = {
                        let mut p = vec![false; cutoff_len];
                        for &i in &a_nonzero {
                            p[i] = true;
                        }
                        p
                    };
                    let bp = {
                        let mut p = vec![false; cutoff_len];
                        for &i in &b_nonzero {
                            p[i] = true;
                        }
                        p
                    };
                    let expected =
                        reference(&a, a_coefficients, &ap, b_coefficients, &bp);

                    for threads in [1usize, 4usize] {
                        let pool = rayon::ThreadPoolBuilder::new()
                            .num_threads(threads)
                            .build()
                            .expect("pool");
                        pool.install(|| {
                            // Public direct kernel.
                            let mut out =
                                vec![NotNan::<f64>::default(); len];
                            LieSeries::commutator_coefficients_with_nonzero(
                                &a,
                                a_coefficients,
                                &a_nonzero,
                                b_coefficients,
                                &b_nonzero,
                                &mut out,
                            );
                            assert_eq!(
                                out, expected,
                                "public kernel differs (d={d}, m={m}, t={threads}, \
                                 a={:?}, b={:?})",
                                a_coefficients.iter().map(|x| x.into_inner()).collect::<Vec<_>>(),
                                b_coefficients.iter().map(|x| x.into_inner()).collect::<Vec<_>>()
                            );

                            // Class-contiguous kernel: operands class-
                            // ordered, result mapped back through `inv`.
                            let order = a.class_order();
                            let a_cls = a.class_coefficients(&order, a_coefficients);
                            let b_cls = b.class_coefficients(&order, b_coefficients);
                            let a_nz_cls: Vec<usize> =
                                a_nonzero.iter().map(|&i| order.inv()[i] as usize).collect();
                            let b_nz_cls: Vec<usize> =
                                b_nonzero.iter().map(|&i| order.inv()[i] as usize).collect();
                            let mut out_cls =
                                vec![NotNan::<f64>::default(); len];
                            a.class_commutation(
                                &order,
                                &a_cls,
                                &a_nz_cls,
                                &b_cls,
                                &b_nz_cls,
                                &mut out_cls,
                                &mut GatingCache::default(),
                            );
                            let public = a.public_coefficients(&order, &out_cls);
                            assert_eq!(
                                public, expected,
                                "class kernel differs (d={d}, m={m}, t={threads})"
                            );
                        });
                    }
                }
            }
        }
    }

    /// ADVERSARIAL (write-access division): a word class stays OWNED under
    /// total cancellation — `a == b` makes every term exactly zero, yet
    /// every reachable word is still swept (one single-writer store of the
    /// exact zero) and still reported by the collecting variant: the
    /// reported targets are a WRITE-set superset, never a value-filtered
    /// set. Distinct failure: guards target ownership and the superset
    /// property, not arithmetic (the fan-in test above).
    #[test]
    fn write_class_canceled_word_stays_owned_and_reported() {
        use lyndon_rs::lyndon::{LyndonBasis, Sort};
        use ordered_float::NotNan;

        for (d, m) in [(2usize, 5usize), (3usize, 4usize)] {
            let words: Vec<LyndonWord<u8>> =
                LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
            let len = words.len();
            let coeffs: Vec<_> = (0..len)
                .map(|i| NotNan::new(((i % 9 + 1) as f64) / 4.0 - 1.1).unwrap())
                .collect();
            let a: LieSeries<u8, NotNan<f64>> =
                LieSeries::new(words.clone(), coeffs.clone());
            let b: LieSeries<u8, NotNan<f64>> = LieSeries::new(words, coeffs);
            let full: Vec<usize> = (0..len).collect();

            // [A, A] = 0 exactly (float products are commutative), but
            // every word with a feasible decomposition is still written.
            let mut out = vec![NotNan::<f64>::default(); len];
            LieSeries::commutator_coefficients_with_nonzero(
                &a,
                &a.coefficients,
                &full,
                &b.coefficients,
                &full,
                &mut out,
            );
            for (w, v) in out.iter().enumerate() {
                assert_eq!(*v, NotNan::<f64>::default(), "word {w} not exactly zero");
            }

            // The collecting variant must still report every written word:
            // value-filtered targets would be unsound for list reuse (a
            // canceled target can go nonzero again when values change).
            let mut dirty = vec![0u64; len / 64 + 1];
            let mut targets = Vec::new();
            LieSeries::commutator_coefficients_with_nonzero_collecting(
                &a,
                &a.coefficients,
                &full,
                &b.coefficients,
                &full,
                &mut out,
                &mut dirty,
                &mut targets,
            );
            // Independent expectation: every word targeted by any active
            // contribution of any entry (full supports -> every entry).
            let table = &a.feasible_decompositions;
            let entries = table.entries();
            let decomp = table.decomp_indices_rel();
            let degrees = table.index_degrees();
            let mut expected: Vec<usize> = Vec::new();
            for (ei, entry) in entries[..entries.len() - 1].iter().enumerate() {
                let (i, j) = (entry.i as usize, entry.j as usize);
                let from = entry.decomp_start as usize;
                let to = entries[ei + 1].decomp_start as usize;
                let rs = table.degree_start(degrees[i] as usize + degrees[j] as usize);
                for &rel in &decomp[from..to] {
                    expected.push(rs + rel as usize);
                }
            }
            expected.sort_unstable();
            expected.dedup();
            // The sink reports kernel-VISIBLE positions only: the degree-
            // `max_degree` tail lies above the support cutoff (nothing can
            // consume those values in a truncated BCH fold), so both sides
            // filter at `degree_start(max_degree)`.
            let cutoff = table.degree_start(table.max_degree());
            expected.retain(|&w| w < cutoff);
            targets.sort_unstable();
            assert_eq!(
                targets, expected,
                "collecting targets lost canceled words (write-set superset violated)"
            );
        }
    }

    /// Bundling decides only which thread runs which anagram unit; the
    /// per-word accumulation order is unchanged, so the result must be
    /// bit-identical for any thread count. Guards the work-balanced bundle
    /// builder: the parallel sweep must never reorder, duplicate, or drop
    /// an accumulation.
    #[test]
    fn commutator_is_bit_identical_across_thread_counts() {
        use lyndon_rs::lyndon::{LyndonBasis, Sort};
        use ordered_float::NotNan;

        for (d, m) in [(2usize, 12usize), (3, 8), (4, 6)] {
            let words: Vec<LyndonWord<u8>> =
                LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
            let a_coefficients: Vec<_> = (0..words.len())
                .map(|i: usize| NotNan::new(((i % 11 + 1) as f64) / 7.0 - 0.9).unwrap())
                .collect();
            let b_coefficients: Vec<_> = (0..words.len())
                .map(|i: usize| NotNan::new(((i * 5 + 3) % 13) as f64 / 6.0 - 1.3).unwrap())
                .collect();
            let a: LieSeries<u8, NotNan<f64>> = LieSeries::new(words.clone(), a_coefficients);
            let b: LieSeries<u8, NotNan<f64>> = LieSeries::new(words, b_coefficients);

            let serial = {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(1)
                    .build()
                    .expect("1-thread pool");
                let mut out = vec![NotNan::<f64>::default(); a.coefficients.len()];
                pool.install(|| {
                    LieSeries::commutator_coefficients(
                        &a,
                        &a.coefficients,
                        &b.coefficients,
                        &mut out,
                    )
                });
                out
            };
            for threads in [2usize, 4usize, 8usize] {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .expect("thread pool");
                let mut out = vec![NotNan::<f64>::default(); a.coefficients.len()];
                pool.install(|| {
                    LieSeries::commutator_coefficients(
                        &a,
                        &a.coefficients,
                        &b.coefficients,
                        &mut out,
                    )
                });
                assert_eq!(
                    serial, out,
                    "parallel result differs from serial for d={d}, m={m}, threads={threads}"
                );
            }
        }
    }

    /// The trait's internal working mode (`class_commutation`) must agree
    /// bit-for-bit with the direct kernel after one
    /// `public_coefficients` epilogue, at every thread count; the
    /// collecting variant must report the class-indexed image of the
    /// direct layout's first-touch sequence.
    #[test]
    fn class_commutation_round_trip_is_bit_identical() {
        use lyndon_rs::lyndon::{LyndonBasis, Sort};
        use ordered_float::NotNan;

        use super::ClassOrderedCommutation;

        for (d, m, force) in [(2usize, 4usize, true), (3usize, 10usize, false)] {
            let words: Vec<LyndonWord<u8>> =
                LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
            let a_coefficients: Vec<_> = (0..words.len())
                .map(|i: usize| NotNan::new(((i % 11 + 1) as f64) / 7.0 - 0.9).unwrap())
                .collect();
            let b_coefficients: Vec<_> = (0..words.len())
                .map(|i: usize| NotNan::new(((i * 5 + 3) % 13) as f64 / 6.0 - 1.3).unwrap())
                .collect();
            let mut direct = LieSeries::new(words.clone(), a_coefficients.clone());
            let class = LieSeries::new(words, b_coefficients.clone());
            if force {
                Arc::make_mut(&mut direct.feasible_decompositions).clear_class_order();
            }
            let order = class.class_order();
            assert_eq!(order.inv().len(), class.coefficients.len());

            let a_cls = class.class_coefficients(&order, &a_coefficients);
            let b_cls = class.class_coefficients(&order, &b_coefficients);
            let inv_of = |i: usize| order.inv()[i] as usize;
            let a_nonzero = direct.nonzero_coefficient_indices(&a_coefficients);
            let b_nonzero = direct.nonzero_coefficient_indices(&b_coefficients);
            let a_nz_cls: Vec<usize> = a_nonzero.iter().copied().map(inv_of).collect();
            let b_nz_cls: Vec<usize> = b_nonzero.iter().copied().map(inv_of).collect();

            let run = |threads: usize, series: &LieSeries<u8, NotNan<f64>>| {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .expect("thread pool");
                let mut out = vec![NotNan::<f64>::default(); series.coefficients.len()];
                pool.install(|| {
                    LieSeries::commutator_coefficients(
                        series,
                        &a_coefficients,
                        &b_coefficients,
                        &mut out,
                    )
                });
                out
            };

            let reference = run(1, &direct);
            for threads in [1usize, 2, 4, 8] {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .expect("thread pool");
                let mut result = vec![NotNan::<f64>::default(); class.coefficients.len()];
                pool.install(|| {
                    class.class_commutation(
                        &order,
                        &a_cls,
                        &a_nz_cls,
                        &b_cls,
                        &b_nz_cls,
                        &mut result,
                        &mut GatingCache::default(),
                    )
                });
                let public = class.public_coefficients(&order, &result);
                assert_eq!(
                    reference, public,
                    "class round trip differs for d={d}, m={m}, threads={threads}"
                );
            }

            // Collecting variant: class-indexed first-touch sequence.
            let mut direct_result = vec![NotNan::<f64>::default(); direct.coefficients.len()];
            let mut direct_dirty = vec![0u64; direct.coefficients.len() / 64 + 1];
            let mut direct_targets = Vec::new();
            LieSeries::commutator_coefficients_with_nonzero_collecting(
                &direct,
                &a_coefficients,
                &a_nonzero,
                &b_coefficients,
                &b_nonzero,
                &mut direct_result,
                &mut direct_dirty,
                &mut direct_targets,
            );
            let mut class_result = vec![NotNan::<f64>::default(); class.coefficients.len()];
            let mut class_dirty = vec![0u64; class.coefficients.len() / 64 + 1];
            let mut class_targets = Vec::new();
            class.class_commutation_with_nonzero_collecting(
                &order,
                &a_cls,
                &a_nz_cls,
                &b_cls,
                &b_nz_cls,
                &mut class_result,
                &mut class_dirty,
                &mut class_targets,
            );
            for (k, &src) in order.inv().iter().enumerate() {
                assert_eq!(
                    direct_result[k], class_result[src as usize],
                    "collecting result differs for d={d}, m={m} at public {k}"
                );
            }
            // Class targets are internal positions: map them back with
            // `perm` (internal -> public) before comparing. Under the
            // write-access division each variant emits its targets sorted
            // ascending in ITS OWN layout's position order (one store per
            // word class, classes in position order), so the public
            // sequences differ by the layout permutation — the invariant is
            // the same target SET, each position reported exactly once
            // (per-word accumulation order is guarded by the bit-identical
            // result comparison above).
            let perm = order.perm();
            let mut relabeled: Vec<usize> = class_targets
                .iter()
                .copied()
                .map(|p| perm[p] as usize)
                .collect();
            relabeled.sort_unstable();
            let mut direct_sorted = direct_targets.clone();
            direct_sorted.sort_unstable();
            assert_eq!(direct_sorted, relabeled, "collecting targets differ");
        }
    }

    /// The class-contiguous scatter layout must produce bit-identical
    /// results to the direct layout, at every thread count, on both kernel
    /// entry points: the layout is a pure relabeling of write addresses and
    /// never reorders the per-word accumulation sequence. `(2, 4)` forces
    /// the layout on a small table (fast build); `(3, 10)` qualifies
    /// through the real slice-size threshold (its degree-10 slice exceeds
    /// 4096 words).
    #[test]
    fn commutator_is_bit_identical_across_scatter_layouts() {
        use lyndon_rs::lyndon::{LyndonBasis, Sort};
        use ordered_float::NotNan;

        for (d, m, force) in [(2usize, 4usize, true), (3usize, 10usize, false)] {
            let words: Vec<LyndonWord<u8>> =
                LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
            let a_coefficients: Vec<_> = (0..words.len())
                .map(|i: usize| NotNan::new(((i % 11 + 1) as f64) / 7.0 - 0.9).unwrap())
                .collect();
            let b_coefficients: Vec<_> = (0..words.len())
                .map(|i: usize| NotNan::new(((i * 5 + 3) % 13) as f64 / 6.0 - 1.3).unwrap())
                .collect();
            let mut direct = LieSeries::new(words.clone(), a_coefficients.clone());
            let mut class = LieSeries::new(words, b_coefficients.clone());
            if force {
                Arc::make_mut(&mut class.feasible_decompositions).force_class_order();
            } else {
                // (3, 10) auto-qualifies through the threshold: strip the
                // reference's layout so the pair compares direct vs class
                // on the same table.
                Arc::make_mut(&mut direct.feasible_decompositions).clear_class_order();
            }
            assert!(
                direct
                    .feasible_decompositions
                    .cached_class_order()
                    .is_none(),
                "expected the reference series to run the direct layout"
            );
            assert!(
                class.feasible_decompositions.cached_class_order().is_some(),
                "expected the layout series to carry a class order"
            );
            let len = direct.coefficients.len();

            let run = |series: &LieSeries<u8, NotNan<f64>>, threads: usize| {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .expect("thread pool");
                let mut out = vec![NotNan::<f64>::default(); len];
                pool.install(|| {
                    LieSeries::commutator_coefficients(
                        series,
                        &a_coefficients,
                        &b_coefficients,
                        &mut out,
                    )
                });
                out
            };

            let reference = run(&direct, 1);
            for threads in [1usize, 2, 4, 8] {
                let out = run(&class, threads);
                assert_eq!(
                    reference, out,
                    "class layout differs from direct for d={d}, m={m}, threads={threads}"
                );
            }

            // Collecting entry point: result, nonzero set, and first-touch
            // order must all match the direct layout exactly.
            let collect = |series: &LieSeries<u8, NotNan<f64>>| {
                let a_nonzero = series.nonzero_coefficient_indices(&a_coefficients);
                let b_nonzero = series.nonzero_coefficient_indices(&b_coefficients);
                let mut result = vec![NotNan::<f64>::default(); len];
                let mut dirty = vec![0u64; len / 64 + 1];
                let mut targets = Vec::new();
                LieSeries::commutator_coefficients_with_nonzero_collecting(
                    series,
                    &a_coefficients,
                    &a_nonzero,
                    &b_coefficients,
                    &b_nonzero,
                    &mut result,
                    &mut dirty,
                    &mut targets,
                );
                (result, targets)
            };
            let (direct_result, direct_targets) = collect(&direct);
            let (class_result, class_targets) = collect(&class);
            assert_eq!(direct_result, class_result, "collecting result differs");
            assert_eq!(direct_targets, class_targets, "collecting targets differ");
        }
    }
}

#[cfg(test)]
mod anagram {
    use super::*;
    use lyndon_rs::lyndon::{LyndonBasis, Sort};
    use ordered_float::NotNan;

    /// The free Lie algebra is multigraded by letter content: the bracket of
    /// two content-homogeneous elements is content-homogeneous of the summed
    /// multiset, and the Lyndon basis refines the grading. So every
    /// decomposition word of `[basis[i], basis[j]]` must carry exactly the
    /// letter multiset of `basis[i]` + `basis[j]` — the invariant that lets
    /// the commutation kernel parallelize over target anagram classes with
    /// disjoint writes.
    #[test]
    fn decompositions_are_content_homogeneous() {
        for (d, m) in [(2usize, 8usize), (3, 8), (4, 8)] {
            let words = LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
            let mut letters: Vec<u8> = words
                .iter()
                .filter(|w| w.len() == 1)
                .map(|w| w.letters[0])
                .collect();
            letters.sort_unstable();
            let content = |w: &LyndonWord<u8>| -> Vec<u8> {
                let mut c = vec![0u8; letters.len()];
                for l in &w.letters {
                    let k = letters.iter().position(|a| a == l).unwrap();
                    c[k] += 1;
                }
                c
            };
            let contents: Vec<Vec<u8>> = words.iter().map(content).collect();

            let series = LieSeries::<u8, NotNan<f64>>::new(words, Vec::new());
            let d_len = series.commutator_basis.len();
            for i in 0..d_len {
                for j in (i + 1)..d_len {
                    if let Some((idx, _, _)) = series.decomposition(i, j) {
                        let mut target = contents[i].clone();
                        for k in 0..letters.len() {
                            target[k] += contents[j][k];
                        }
                        for &x in idx {
                            assert_eq!(
                                contents[x as usize], target,
                                "d={d} m={m}: pair ({i}, {j}) decomposition word {x}                                  violates content homogeneity"
                            );
                        }
                    }
                }
            }
        }
    }

    /// The raw-float fast path (`NotNan<f64>` / `f64` dispatch) and the
    /// generic path (`Ratio<i128>`) must agree. Integer-valued coefficients
    /// keep every intermediate exactly representable in `f64` (magnitudes
    /// stay far below 2^53), so the comparison is exact. (Ported from the
    /// original raw-float fast path change.)
    #[test]
    fn raw_float_kernel_matches_rationals() {
        use lyndon_rs::lyndon::{LyndonBasis, Sort};
        use num_rational::Ratio;
        use num_traits::ToPrimitive;

        for (d, m) in [(2usize, 6usize), (3, 5)] {
            let words = LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
            let coeffs = |salt: usize| {
                (0..words.len())
                    .map(|i| ((i * 7 + salt * 13) % 21) as i128 - 10)
                    .collect::<Vec<_>>()
            };
            let (a_int, b_int) = (coeffs(1), coeffs(2));

            // Raw-float path.
            let a_f = LieSeries::<u8, NotNan<f64>>::new(
                words.clone(),
                a_int
                    .iter()
                    .map(|&x| NotNan::new(x as f64).unwrap())
                    .collect::<Vec<_>>(),
            );
            let b_f = LieSeries::<u8, NotNan<f64>>::new(
                words.clone(),
                b_int
                    .iter()
                    .map(|&x| NotNan::new(x as f64).unwrap())
                    .collect::<Vec<_>>(),
            );
            let ab_f = a_f.commutator(&b_f);

            // Generic exact path.
            let a_r = LieSeries::<u8, Ratio<i128>>::new(
                words.clone(),
                a_int.iter().map(|&x| Ratio::from_integer(x)).collect(),
            );
            let b_r = LieSeries::<u8, Ratio<i128>>::new(
                words.clone(),
                b_int.iter().map(|&x| Ratio::from_integer(x)).collect(),
            );
            let ab_r = a_r.commutator(&b_r);

            for (x, y) in ab_f.coefficients.iter().zip(&ab_r.coefficients) {
                assert_eq!(
                    x.into_inner(),
                    y.to_f64().unwrap(),
                    "d={d} m={m}: raw float vs exact rationals"
                );
            }
        }
    }

    /// The raw helpers must behave bitwise like the primitive float ops and
    /// must NOT panic where the wrapper's arithmetic would: overflow to
    /// infinity and its NaN cancellation are the caller's audit (NaN policy
    /// of `raw_mul`).
    #[test]
    fn raw_ops_match_primitive_semantics_without_panic() {
        use num_rational::Ratio;

        // Overflow: the wrapper's Mul panics on NaN results only; plain
        // overflow to inf is fine in both. Check bitwise equality with the
        // primitive for representative inputs, including the NaN-producing
        // combination the checked path would panic on.
        let cases = [
            (3.0f64, 5.0f64),
            (-2.5, 4.25),
            (f64::MAX, f64::MAX),   // -> inf
            (f64::MAX, -f64::MAX),  // -> -inf
        ];
        for (x, y) in cases {
            let a = NotNan::new(x).unwrap();
            let b = NotNan::new(y).unwrap();
            let raw = raw_mul(&a, &b);
            assert_eq!(raw.into_inner().to_bits(), (x * y).to_bits());
        }
        // NaN cancellation: inf + (-inf) — the wrapper's AddAssign panics,
        // the raw helper writes the NaN (audit is the caller's job).
        let mut acc = NotNan::new(f64::INFINITY).unwrap();
        let neg_inf = NotNan::new(f64::NEG_INFINITY).unwrap();
        raw_add_assign(&mut acc, &neg_inf);
        assert!(is_nan_f64(acc.into_inner()));

        // f32 mirrors f64.
        let a = NotNan::new(f32::MAX).unwrap();
        let raw = raw_mul(&a, &a);
        assert_eq!(raw.into_inner().to_bits(), (f32::MAX * f32::MAX).to_bits());

        // The generic (non-float) path is untouched: exact multiplication.
        let r = Ratio::new(7i128, 3);
        let s = Ratio::new(-2i128, 5);
        let mut acc_r = r.clone();
        raw_add_assign(&mut acc_r, &s);
        assert_eq!(acc_r, r + s);
        assert_eq!(raw_mul(&r, &s), r * s);
    }

    #[inline]
    fn is_nan_f64(x: f64) -> bool {
        // `x != x` trips clippy::eq_op (deny-by-default); is_nan is the
        // same predicate for f64.
        x.is_nan()
    }
}

#[cfg(test)]
mod obj1_probe {
    use super::*;
    use lyndon_rs::lyndon::{LyndonBasis, Sort};
    use ordered_float::NotNan;

    #[test]
    #[ignore = "probe"]
    fn probe_entry_contribution_ratio() {
        for (d, m) in [(2usize, 12usize), (3usize, 8usize), (4usize, 8usize), (2usize, 8usize), (3usize, 12usize)] {
            let words: Vec<LyndonWord<u8>> =
                LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
            let len = words.len();
            let a: LieSeries<u8, NotNan<f64>> =
                LieSeries::new(words, vec![NotNan::new(1.0).unwrap(); len]);
            let table = &a.feasible_decompositions;
            let entries = table.entries();
            let e = entries.len() - 1;
            let c: usize = (0..e)
                .map(|ei| {
                    entries[ei + 1].decomp_start as usize - entries[ei].decomp_start as usize
                })
                .sum();
            let rows_gt1 = (0..e)
                .filter(|&ei| {
                    entries[ei + 1].decomp_start - entries[ei].decomp_start > 1
                })
                .count();
            println!(
                "d={d} m={m}: basis={len} entries={e} contribs={c} rho={:.3} rows_gt1={rows_gt1} ({:.1}%)",
                c as f64 / e as f64,
                100 * rows_gt1 / e.max(1)
            );
        }
    }
}
