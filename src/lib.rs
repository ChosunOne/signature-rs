use std::fmt::Debug;
use std::hash::Hash;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, RemAssign, Sub, SubAssign};

use num_traits::{FromPrimitive, Num};

use crate::generators::Generator;
use crate::lie_series::LieSeries;

pub mod bch;
pub mod bch_series_generator;
pub mod commutator;
pub mod constants;
pub mod generators;
pub mod lie_series;
pub mod log_sig;
pub mod lyndon;
pub mod rooted_tree;

pub trait Arith:
    Clone
    + Default
    + Debug
    + Sized
    + Add<Output = Self>
    + AddAssign
    + Sub<Output = Self>
    + SubAssign
    + Mul<Output = Self>
    + MulAssign
    + Div<Output = Self>
    + DivAssign
    + RemAssign
    + Neg<Output = Self>
    + PartialOrd
    + PartialEq
    + Eq
    + Ord
    + Sized
    + Num
    + FromPrimitive
    + Hash
{
}

impl<
    T: Clone
        + Default
        + Debug
        + Sized
        + Add<Output = Self>
        + AddAssign
        + Sub<Output = Self>
        + SubAssign
        + Mul<Output = Self>
        + MulAssign
        + Div<Output = Self>
        + DivAssign
        + RemAssign
        + Neg<Output = Self>
        + PartialOrd
        + Ord
        + PartialEq
        + Eq
        + Sized
        + Num
        + FromPrimitive
        + Hash,
> Arith for T
{
}

pub trait BCHCoefficientGenerator {
    fn generate_bch_coefficients<U: Arith + Send + Sync>(&self) -> Vec<U>;
}

pub trait LieSeriesGenerator<T: Generator, U: Arith + Send + Sync> {
    fn generate_lie_series(&self) -> LieSeries<T, U>;
}

fn binomial<U: Arith>(n: usize, k: usize) -> U {
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

#[must_use]
pub fn bernoulli_sequence<U: Arith>(max_n: usize) -> Vec<U> {
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
