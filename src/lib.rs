use std::fmt::Debug;
use std::hash::Hash;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, RemAssign, Sub, SubAssign};

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::Ratio;

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
    + Ord
    + PartialEq
    + Eq
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
        + Eq
        + PartialEq
        + Hash,
> Arith for T
{
}

pub trait Int: Arith + Integer + From<i8> + From<i32> + From<i64> + From<i128> + From<u64> {}

impl Int for i128 {}
impl Int for BigInt {}

pub trait BCHCoefficientGenerator {
    fn generate_bch_coefficients<U: Int + Send + Sync>(&self) -> Vec<Ratio<U>>;
}

pub trait LieSeriesGenerator<
    const N: usize,
    T: Generator<Letter = T>,
    U: Int + Hash + Arith + Send + Sync,
>
{
    fn generate_lie_series(&self) -> LieSeries<N, T, U>;
}

fn binomial<U: Int>(n: usize, k: usize) -> U {
    if k == 0 {
        return U::from(1);
    }

    let k = k.min(n - k);
    let mut num = U::from(1);
    let mut den = U::from(1);

    for i in 0..k {
        num *= U::from((n - i) as u64);
        den *= U::from(i as u64 + 1);
    }

    num / den
}

#[must_use]
pub fn bernoulli_sequence<U: Int>(max_n: usize) -> Vec<Ratio<U>> {
    let mut b = vec![Ratio::default(); max_n + 1];

    b[0] = Ratio::new(1.into(), 1.into());
    if max_n > 0 {
        b[1] = Ratio::new((-1).into(), 2.into());
    }

    for n in 2..=max_n {
        if n % 2 == 1 {
            continue;
        }
        let mut sum = Ratio::default();
        for k in 0..n {
            sum += Ratio::from(binomial::<U>(n + 1, k)) * &b[k];
        }

        b[n] = -sum * Ratio::new(U::from(1), (n as u64 + 1).into());
    }

    // Use positive sign convention
    if max_n > 0 {
        b[1] = Ratio::new(1.into(), 2.into());
    }

    b
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_bernoulli_sequence() {
        let seq = bernoulli_sequence::<i128>(10);
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
