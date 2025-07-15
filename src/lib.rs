use std::ops::{AddAssign, DivAssign, MulAssign, Neg, RemAssign, SubAssign};

#[cfg(feature = "progress")]
use indicatif::ProgressBar;

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::Ratio;

pub mod bch;
pub mod constants;
pub mod lie_series;
pub mod lyndon;
pub mod rooted_tree;

pub trait Int:
    Clone
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
    fn generate_bch_coefficients<U: Int>(&self) -> Vec<Ratio<U>>;
}

fn binomial<
    U: Clone
        + Integer
        + AddAssign
        + SubAssign
        + MulAssign
        + RemAssign
        + DivAssign
        + From<i32>
        + From<u64>,
>(
    n: usize,
    k: usize,
) -> Ratio<U> {
    if k == 0 {
        return Ratio::new(1.into(), 1.into());
    }

    let k = k.min(n - k);
    let mut result = Ratio::new(1.into(), 1.into());

    for i in 0..k {
        result *= Ratio::new(((n - i) as u64).into(), (i as u64 + 1).into());
    }

    result
}

pub fn bernoulli_sequence<
    U: Clone
        + Integer
        + AddAssign
        + SubAssign
        + MulAssign
        + RemAssign
        + DivAssign
        + Neg<Output = U>
        + From<i32>
        + From<u64>,
>(
    max_n: usize,
) -> Vec<Ratio<U>> {
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
            sum += binomial(n + 1, k) * &b[k];
        }

        b[n] = -sum * Ratio::new(1.into(), (n as u64 + 1).into());
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
