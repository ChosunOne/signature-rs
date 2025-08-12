use num_rational::Ratio;

use crate::{
    Arith, Int, LieSeriesGenerator,
    bch_series_generator::BchSeriesGenerator,
    comm,
    commutator::{Commutator, CommutatorTerm},
    lie_series::LieSeries,
    lyndon::{Generator, Lexicographical, LyndonBasis, LyndonWord},
};
use std::{
    collections::HashMap,
    hash::Hash,
    ops::{Index, IndexMut},
};

#[derive(Debug, Clone)]
pub struct LogSignature<
    const N: usize,
    T: Generator<Letter = T> + Send + Sync,
    U: Int + Arith + Hash + Send + Sync,
> {
    series: LieSeries<N, T, U>,
    bch_series: LieSeries<2, u8, U>,
}

impl<const N: usize, T: Generator<Letter = T> + Send + Sync, U: Int + Hash + Send + Sync>
    Index<usize> for LogSignature<N, T, U>
{
    type Output = Ratio<U>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.series[index]
    }
}

impl<const N: usize, T: Generator<Letter = T> + Send + Sync, U: Int + Hash + Send + Sync>
    Index<LyndonWord<N, T>> for LogSignature<N, T, U>
{
    type Output = Ratio<U>;

    fn index(&self, index: LyndonWord<N, T>) -> &Self::Output {
        &self.series[index]
    }
}

impl<const N: usize, T: Generator<Letter = T> + Send + Sync, U: Int + Hash + Send + Sync>
    Index<&LyndonWord<N, T>> for LogSignature<N, T, U>
{
    type Output = Ratio<U>;

    fn index(&self, index: &LyndonWord<N, T>) -> &Self::Output {
        &self.series[index]
    }
}

impl<const N: usize, T: Generator<Letter = T> + Send + Sync, U: Int + Hash + Send + Sync>
    IndexMut<usize> for LogSignature<N, T, U>
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.series[index]
    }
}

impl<const N: usize, T: Generator<Letter = T> + Send + Sync, U: Int + Hash + Send + Sync>
    IndexMut<LyndonWord<N, T>> for LogSignature<N, T, U>
{
    fn index_mut(&mut self, index: LyndonWord<N, T>) -> &mut Self::Output {
        &mut self.series[index]
    }
}

impl<const N: usize, T: Generator<Letter = T> + Send + Sync, U: Int + Hash + Send + Sync>
    IndexMut<&LyndonWord<N, T>> for LogSignature<N, T, U>
{
    fn index_mut(&mut self, index: &LyndonWord<N, T>) -> &mut Self::Output {
        &mut self.series[index]
    }
}

impl<const N: usize, T: Generator<Letter = T> + Send + Sync, U: Int + Hash + Send + Sync>
    LogSignature<N, T, U>
{
    #[must_use]
    pub fn new(max_lyndon_word_length: usize) -> Self {
        let basis = LyndonBasis::<N, T, Lexicographical>::generate_basis(max_lyndon_word_length);
        let bch_basis =
            LyndonBasis::<2, u8, Lexicographical>::generate_basis(max_lyndon_word_length);
        let bch_series_generator = BchSeriesGenerator::<2, u8>::new(bch_basis.clone());
        let bch_series = bch_series_generator.generate_lie_series();
        let coefficients = vec![Ratio::<U>::from(U::from(0)); basis.len()];
        let series = LieSeries::<N, T, U>::new(basis, coefficients);

        Self { series, bch_series }
    }

    #[must_use]
    pub fn concatenate(&self, rhs: &Self) -> Self {
        let mut computed_commutations = HashMap::new();
        let mut concatenated_log_sig = self.clone();

        for (i, term) in self.bch_series.commutator_basis.iter().enumerate().skip(1) {
            concatenated_log_sig.series += evaluate_commutator_term(
                term,
                &[&self.series, &rhs.series],
                &mut computed_commutations,
            ) * self.bch_series[i].clone();
        }

        concatenated_log_sig
    }
}

fn evaluate_commutator_term<
    const N: usize,
    T: Generator<Letter = T>,
    U: Int + Hash + Arith + Send + Sync,
>(
    term: &CommutatorTerm<U, u8>,
    series: &[&LieSeries<N, T, U>],
    computed_commutations: &mut HashMap<CommutatorTerm<U, u8>, LieSeries<N, T, U>>,
) -> LieSeries<N, T, U> {
    match term {
        &CommutatorTerm::Atom(a) => series[a as usize].clone(),
        t @ CommutatorTerm::Expression(e) => {
            if computed_commutations.contains_key(t) {
                return computed_commutations[t].clone();
            }
            let a = evaluate_commutator_term(&e.left, series, computed_commutations);
            let b = evaluate_commutator_term(&e.right, series, computed_commutations);
            let result = comm![a, b];
            computed_commutations.insert(t.clone(), result.clone());
            result
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_log_sig_concat() {
        let mut a = LogSignature::<2, u8, i128>::new(3);
        let mut b = LogSignature::new(3);
        a.series.coefficients = [1, 2, 3, 4, 5].map(Ratio::from_integer).to_vec();
        b.series.coefficients = [6, 7, 8, 9, 10].map(Ratio::from_integer).to_vec();
        let c = a.concatenate(&b);
        let expected_coefficients = [
            Ratio::new(7, 1),
            Ratio::new(9, 1),
            Ratio::new(17, 2),
            Ratio::new(121, 12),
            Ratio::new(185, 12),
        ];
        assert_eq!(c.series.coefficients, expected_coefficients);
    }
}
