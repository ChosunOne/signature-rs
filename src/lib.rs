use std::fmt::Debug;
use std::hash::Hash;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, RemAssign, Sub, SubAssign};

use num_traits::{FromPrimitive, Num};

use crate::lie_series::LieSeries;
use crate::lyndon::Generator;

pub mod bch;
pub mod bch_series_generator;
pub mod commutator;
pub mod constants;
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
    + Hash
    + Sized
    + Num
    + FromPrimitive
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
        + Eq
        + PartialEq
        + Hash
        + Sized
        + Num
        + FromPrimitive,
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
        b[1] = -U::one() / U::from_u8(2).expect("Failed to convert from i8");
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

    use super::*;

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
