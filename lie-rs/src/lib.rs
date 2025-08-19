mod bch;
pub mod bch_series_generator;
mod constants;
pub mod lie_series;
pub mod rooted_tree;

use std::ops::{AddAssign, Div, MulAssign, Neg};

pub use lie_series::LieSeries;
use num_traits::{FromPrimitive, One, Zero};

/// Trait for types that can generate Lie series.
///
/// A Lie series is a formal power series in non-commuting variables that
/// respects the Lie algebra structure. This trait provides a uniform
/// interface for generating such series from various sources.
pub trait LieSeriesGenerator<T, U> {
    /// Generates a Lie series from this generator.
    fn generate_lie_series(&self) -> LieSeries<T, U>;
}

/// Computes the binomial coefficient "n choose k" for generic numeric types.
///
/// This function efficiently calculates C(n,k) = n!/(k!(n-k)!) using
/// the multiplicative formula to avoid computing large factorials directly.
pub(crate) fn binomial<U: One + Zero + FromPrimitive + MulAssign + Div<Output = U>>(
    n: usize,
    k: usize,
) -> U {
    if k == 0 {
        return U::one();
    }

    if k > n {
        return U::zero();
    }

    let k = k.min(n - k);
    let mut num = U::one();
    let mut den = U::one();

    for i in 0..k {
        num *= U::from_usize(n - i).expect("Failed to convert from usize");
        den *= U::from_usize(i + 1).expect("Failed to convert from usize");
    }

    num / den
}

/// Generates a sequence of Bernoulli numbers up to the specified maximum index.
///
/// Bernoulli numbers B_n are rational numbers that appear in the series expansions
/// of many important mathematical functions. This implementation uses the recursive
/// formula based on binomial coefficients and applies the positive sign convention
/// for B_1 = +1/2.
#[must_use]
pub(crate) fn bernoulli_sequence<
    U: Clone
        + Default
        + One
        + Zero
        + FromPrimitive
        + Neg<Output = U>
        + Div<Output = U>
        + AddAssign
        + MulAssign,
>(
    max_n: usize,
) -> Vec<U> {
    let mut b = vec![U::default(); max_n + 1];

    b[0] = U::one();
    if max_n > 0 {
        b[1] = -U::one() / U::from_u8(2).expect("Failed to convert from u8");
    }

    for n in 2..=max_n {
        if n % 2 == 1 {
            continue;
        }
        let mut sum = U::default();
        for (k, b_k) in b.iter().enumerate().take(n) {
            sum += binomial::<U>(n + 1, k) * b_k.clone();
        }

        b[n] = -sum / U::from_usize(n + 1).expect("Failed to convert from usize");
    }

    // Use positive sign convention
    if max_n > 0 {
        b[1] = -b[1].clone();
    }

    b
}

#[cfg(test)]
mod test {
    use num_rational::Ratio;
    use ordered_float::NotNan;
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case(0, 0, 1)]
    #[case(5, 0, 1)]
    #[case(100, 0, 1)]
    #[case(1, 1, 1)]
    #[case(2, 2, 1)]
    #[case(3, 3, 1)]
    #[case(4, 4, 1)]
    #[case(5, 5, 1)]
    #[case(5, 6, 0)]
    #[case(4, 1, 4)]
    #[case(4, 2, 6)]
    #[case(4, 3, 4)]
    #[case(4, 4, 1)]
    #[case(30, 15, 155_117_520)]
    fn test_binomial_integer(#[case] n: usize, #[case] k: usize, #[case] expected_result: i128) {
        assert_eq!(binomial::<i128>(n, k), expected_result);
    }

    #[rstest]
    #[case(0, 0, NotNan::from(1))]
    #[case(5, 0, NotNan::from(1))]
    #[case(100, 0, NotNan::from(1))]
    #[case(1, 1, NotNan::from(1))]
    #[case(2, 2, NotNan::from(1))]
    #[case(3, 3, NotNan::from(1))]
    #[case(4, 4, NotNan::from(1))]
    #[case(5, 5, NotNan::from(1))]
    #[case(4, 1, NotNan::from(4))]
    #[case(4, 2, NotNan::from(6))]
    #[case(4, 3, NotNan::from(4))]
    #[case(4, 4, NotNan::from(1))]
    #[case(30, 15, NotNan::from(155_117_520))]
    fn test_binomial_float(
        #[case] n: usize,
        #[case] k: usize,
        #[case] expected_result: NotNan<f64>,
    ) {
        assert_eq!(binomial::<NotNan<f64>>(n, k), expected_result);
    }

    #[rstest]
    #[case(0, 0, Ratio::from_integer(1))]
    #[case(5, 0, Ratio::from_integer(1))]
    #[case(100, 0, Ratio::from_integer(1))]
    #[case(1, 1, Ratio::from_integer(1))]
    #[case(2, 2, Ratio::from_integer(1))]
    #[case(3, 3, Ratio::from_integer(1))]
    #[case(4, 4, Ratio::from_integer(1))]
    #[case(5, 5, Ratio::from_integer(1))]
    #[case(5, 6, Ratio::from_integer(0))]
    #[case(4, 1, Ratio::from_integer(4))]
    #[case(4, 2, Ratio::from_integer(6))]
    #[case(4, 3, Ratio::from_integer(4))]
    #[case(4, 4, Ratio::from_integer(1))]
    #[case(30, 15, Ratio::from_integer(155_117_520))]
    fn test_binomial_ratio(
        #[case] n: usize,
        #[case] k: usize,
        #[case] expected_result: Ratio<i128>,
    ) {
        assert_eq!(binomial::<Ratio<i128>>(n, k), expected_result);
    }

    #[test]
    fn test_bernoulli_sequence() {
        let seq = bernoulli_sequence::<Ratio<i128>>(10);
        let expected_seq = vec![
            Ratio::new(1.into(), 1.into()),
            Ratio::new(1.into(), 2.into()),
            Ratio::new(1.into(), 6.into()),
            Ratio::default(),
            Ratio::new((-1).into(), 30.into()),
            Ratio::default(),
            Ratio::new(1.into(), 42.into()),
            Ratio::default(),
            Ratio::new((-1).into(), 30.into()),
            Ratio::default(),
            Ratio::new(5.into(), 66.into()),
        ];
        assert_eq!(seq.len(), expected_seq.len());
        for (term, expected_term) in seq.iter().zip(&expected_seq) {
            assert_eq!(term, expected_term);
        }
    }
}
