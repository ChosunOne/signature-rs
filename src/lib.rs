use std::fmt::Debug;
use std::ops::{AddAssign, DivAssign, MulAssign, Neg, RemAssign, SubAssign};

#[cfg(feature = "progress")]
use indicatif::ProgressBar;

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::Ratio;

use crate::lie_series::LieSeries;
use crate::lyndon::Generator;

pub mod bch;
pub mod bch_series_generator;
pub mod constants;
pub mod lie_series;
pub mod lyndon;
pub mod rooted_tree;

pub trait Int:
    Clone
    + Default
    + Debug
    + Integer
    + AddAssign
    + SubAssign
    + MulAssign
    + DivAssign
    + RemAssign
    + Neg<Output = Self>
    + From<i8>
    + From<i32>
    + From<i64>
    + From<i128>
    + From<u64>
{
}

impl Int for i128 {}
impl Int for BigInt {}

pub trait BCHCoefficientGenerator {
    fn generate_bch_coefficients<U: Int + Send + Sync>(&self) -> Vec<Ratio<U>>;
}

pub trait LieSeriesGenerator<const N: usize, T: Generator<Letter = T>, U: Int + Send + Sync> {
    fn generate_lie_series(&self) -> LieSeries<N, T, U>;
}

pub trait Commutator<Rhs = Self> {
    type Output;
    /// The commutator operation is represented with `[A, B]` and commonly represents `AB - BA`
    fn commutator(&self, other: Rhs) -> Self::Output;
}

impl<T: Int> Commutator<&Self> for T {
    type Output = T;

    fn commutator(&self, other: &Self) -> Self::Output {
        self.clone() * other.clone() - other.clone() * self.clone()
    }
}

/// Shorthand for applying the commutator operation `[A, B]`.
macro_rules! comm {
    ($a:expr, $b:expr) => {
        $a.commutator(&$b)
    };
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

    #[test]
    fn test_commutators_int() {
        assert_eq!(comm![1, 2], 0);
    }

    #[test]
    fn test_commutators_bigint() {
        let a = BigInt::from(1);
        let b = BigInt::from(2);

        assert_eq!(comm![a, b], BigInt::from(0));
    }
}
