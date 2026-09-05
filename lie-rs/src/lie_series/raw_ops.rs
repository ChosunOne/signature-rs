//! Raw-float fast path for the commutation kernel's coefficient arithmetic.
//!
//! For the float coefficient types used in production the kernel's hot loops
//! run their `*`/`+=` through these helpers instead of the arithmetic impls
//! of the checking [`ordered_float::NotNan`] wrapper (each of which pays a
//! per-operation NaN check, `ucomisd`+`jp` in the sweep's inner loop). The
//! wrapper is `#[repr(transparent)]` over the primitive, so the helpers
//! reinterpret the operands and compute in the primitive type; the
//! `TypeId` comparisons constant-fold per monomorphization, so the dispatch
//! is free for concrete callers and every other coefficient type (e.g.
//! `Ratio<i128>`) keeps the identical generic path.
//!
//! # NaN policy
//! For the raw-float instantiations no per-operation NaN check runs. Finite
//! inputs cannot produce NaN through `*`/`+=`; overflow produces infinities
//! exactly as the checked path does, and only `inf - inf`-style combinations
//! of those can yield NaN. Such a NaN is stored bit-for-bit into the
//! `NotNan` slot (an invalid value for the wrapper's logical invariant, but
//! not undefined behavior — the wrapper has no niche). Callers that persist
//! results in `NotNan` slots audit them once per fold step (the
//! log-signature fold) and fail loudly instead of leaving the broken
//! invariant behind. Results are otherwise bitwise identical to the checked
//! path.
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

#[cfg(test)]
mod tests {
    use super::*;

    /// The raw helpers must behave bitwise like the primitive float ops and
    /// must NOT panic where the wrapper's arithmetic would: overflow to
    /// infinity and its NaN cancellation are the caller's audit (NaN policy
    /// of [`raw_mul`]).
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
            (f64::MAX, f64::MAX),  // -> inf
            (f64::MAX, -f64::MAX), // -> -inf
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
        assert!(acc.into_inner().is_nan());

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
}
