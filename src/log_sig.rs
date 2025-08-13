use num_rational::Ratio;

use crate::{
    Arith, Int, LieSeriesGenerator,
    bch_series_generator::BchSeriesGenerator,
    comm,
    commutator::{Commutator, CommutatorTerm},
    lie_series::LieSeries,
    lyndon::{Generator, LyndonBasis, LyndonWord, Sort},
};
use std::{
    collections::HashMap,
    hash::Hash,
    ops::{Index, IndexMut},
};

#[derive(Debug, Clone)]
pub struct LogSignature<T: Generator + Send + Sync = u8, U: Int + Hash + Send + Sync = i128> {
    series: LieSeries<T, U>,
    bch_series: LieSeries<u8, U>,
}

impl<T: Generator + Send + Sync, U: Int + Hash + Send + Sync> Index<usize> for LogSignature<T, U> {
    type Output = Ratio<U>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.series[index]
    }
}

impl<T: Generator + Send + Sync, U: Int + Hash + Send + Sync> Index<LyndonWord<T>>
    for LogSignature<T, U>
{
    type Output = Ratio<U>;

    fn index(&self, index: LyndonWord<T>) -> &Self::Output {
        &self.series[index]
    }
}

impl<T: Generator + Send + Sync, U: Int + Hash + Send + Sync> Index<&LyndonWord<T>>
    for LogSignature<T, U>
{
    type Output = Ratio<U>;

    fn index(&self, index: &LyndonWord<T>) -> &Self::Output {
        &self.series[index]
    }
}

impl<T: Generator + Send + Sync, U: Int + Hash + Send + Sync> IndexMut<usize>
    for LogSignature<T, U>
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.series[index]
    }
}

impl<T: Generator + Send + Sync, U: Int + Hash + Send + Sync> IndexMut<LyndonWord<T>>
    for LogSignature<T, U>
{
    fn index_mut(&mut self, index: LyndonWord<T>) -> &mut Self::Output {
        &mut self.series[index]
    }
}

impl<T: Generator + Send + Sync, U: Int + Hash + Send + Sync> IndexMut<&LyndonWord<T>>
    for LogSignature<T, U>
{
    fn index_mut(&mut self, index: &LyndonWord<T>) -> &mut Self::Output {
        &mut self.series[index]
    }
}

impl<T: Generator + Send + Sync, U: Int + Hash + Send + Sync> LogSignature<T, U> {
    #[must_use]
    pub fn new(alphabet_size: usize, max_lyndon_word_length: usize) -> Self {
        let basis = LyndonBasis::<T>::new(alphabet_size, Sort::Lexicographical)
            .generate_basis(max_lyndon_word_length);
        let bch_basis = LyndonBasis::new(2, Sort::Lexicographical);
        let bch_series_generator = BchSeriesGenerator::new(bch_basis, max_lyndon_word_length);
        let bch_series = bch_series_generator.generate_lie_series();
        let coefficients = vec![Ratio::<U>::from(U::from(0)); basis.len()];
        let series = LieSeries::<T, U>::new(basis, coefficients);

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

fn evaluate_commutator_term<T: Generator, U: Int + Hash + Send + Sync>(
    term: &CommutatorTerm<U, u8>,
    series: &[&LieSeries<T, U>],
    computed_commutations: &mut HashMap<CommutatorTerm<U, u8>, LieSeries<T, U>>,
) -> LieSeries<T, U> {
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
        let mut a = LogSignature::<u8, i128>::new(2, 3);
        let mut b = LogSignature::new(2, 3);
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
