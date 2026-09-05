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
