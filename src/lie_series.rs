use std::collections::HashMap;
use std::fmt::Debug;
use std::hash::Hash;
use std::ops::{Index, IndexMut};

use num_rational::Ratio;

use crate::{
    Int, comm,
    commutator::{Commutator, CommutatorTerm},
    lyndon::{Generator, LyndonWord},
};

pub struct LieSeries<const N: usize, T: Generator<Letter = T>, U: Hash + Int + Send + Sync> {
    /// The Lyndon basis for the series
    basis: Vec<LyndonWord<N, T>>,
    /// The commutator basis for the series
    commutator_basis: Vec<CommutatorTerm<U, T>>,
    /// A map for converting arbitrary commutator terms to basis elements
    commutator_basis_map: HashMap<CommutatorTerm<U, T>, CommutatorTerm<U, T>>,
    /// A map for locating a given term's index in the basis
    commutator_basis_index_map: HashMap<CommutatorTerm<U, T>, usize>,
    /// The coefficients for each of the terms in the series
    coefficients: Vec<Ratio<U>>,
}

impl<const N: usize, T: Generator<Letter = T>, U: Hash + Int + Send + Sync> Index<usize>
    for LieSeries<N, T, U>
{
    type Output = Ratio<U>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.coefficients[index]
    }
}

impl<const N: usize, T: Generator<Letter = T>, U: Hash + Int + Send + Sync> IndexMut<usize>
    for LieSeries<N, T, U>
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.coefficients[index]
    }
}

impl<const N: usize, T: Generator<Letter = T>, U: Hash + Int + Send + Sync> LieSeries<N, T, U> {
    pub fn new(basis: Vec<LyndonWord<N, T>>, coefficients: Vec<Ratio<U>>) -> Self {
        let mut commutator_basis = Vec::<CommutatorTerm<U, T>>::with_capacity(basis.len());
        for word in basis.iter() {
            commutator_basis.push(CommutatorTerm::from(word));
        }
        let mut commutator_basis_map = HashMap::new();

        let mut commutator_basis_index_map = HashMap::new();
        for (i, a) in commutator_basis.iter().enumerate() {
            commutator_basis_index_map.insert(a.clone(), i);
            match a {
                CommutatorTerm::Expression(e) => {
                    let mut e = e.clone();
                    e.coefficient = -e.coefficient.clone();
                    commutator_basis_index_map.insert(CommutatorTerm::Expression(e), i);
                }
                _ => {}
            }
        }

        for a in commutator_basis.iter() {
            for b in commutator_basis.iter() {
                if a == b {
                    continue;
                }
                let term = comm![a, b];
                let mut lyndon_term = term.clone();
                lyndon_term.lyndon_sort();
                // Only include terms that are in our basis
                if commutator_basis_index_map.contains_key(&lyndon_term) {
                    commutator_basis_map.insert(term, lyndon_term);
                }
            }
        }

        Self {
            basis,
            commutator_basis,
            commutator_basis_map,
            commutator_basis_index_map,
            coefficients,
        }
    }
}

impl<const N: usize, T: Generator<Letter = T> + Debug + Clone, U: Hash + Int + Send + Sync>
    Commutator<&Self> for LieSeries<N, T, U>
{
    type Output = Self;

    /// Calculates the lie bracket `[A, B]` for a lie series for terms within the commutator basis.
    fn commutator(&self, other: &Self) -> Self::Output {
        let mut coefficients = vec![Ratio::<U>::default(); self.coefficients.len()];
        for i in 0..self.coefficients.len() {
            let a = &self.commutator_basis[i];
            for j in 0..other.coefficients.len() {
                if i == j {
                    continue;
                }
                let b = &other.commutator_basis[j];
                let comm_term = comm![a, b];
                let Some(basis_term) = self.commutator_basis_map.get(&comm_term) else {
                    continue;
                };

                let basis_index = self.commutator_basis_index_map[basis_term];
                let CommutatorTerm::Expression(comm_expr) = basis_term else {
                    panic!("Failed to create commutator expression from term");
                };

                coefficients[basis_index] +=
                    self[i].clone() * other[j].clone() * comm_expr.coefficient.clone();
            }
        }
        Self {
            basis: self.basis.clone(),
            commutator_basis: self.commutator_basis.clone(),
            commutator_basis_map: self.commutator_basis_map.clone(),
            commutator_basis_index_map: self.commutator_basis_index_map.clone(),
            coefficients,
        }
    }
}

#[cfg(test)]
mod test {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case(2, vec![1, 2, 3], vec![4, 5, 6], vec![0, 0, -3])]
    #[case(2, vec![3, 2, 1], vec![1, 2, 3], vec![0, 0, 4])]
    #[case(3, vec![1, 2, 3, 4, 5], vec![6, 7, 8, 9, 10], vec![0, 0, -5, -10, 5])]
    fn test_lie_series_commutation(
        #[case] basis_depth: usize,
        #[case] a_coefficients: Vec<i128>,
        #[case] b_coefficients: Vec<i128>,
        #[case] expected_coefficients: Vec<i128>,
    ) {
        use crate::lyndon::{Lexicographical, LyndonBasis};

        let basis = LyndonBasis::<2, u8, Lexicographical>::generate_basis(basis_depth);
        let a_coefficients = a_coefficients
            .into_iter()
            .map(|x| Ratio::<i128>::from_integer(x))
            .collect::<Vec<_>>();
        let a = LieSeries::new(basis.clone(), a_coefficients);
        let b_coefficients = b_coefficients
            .into_iter()
            .map(|x| Ratio::<i128>::from_integer(x))
            .collect::<Vec<_>>();
        let b = LieSeries::new(basis.clone(), b_coefficients);
        let expected_coefficients = expected_coefficients
            .into_iter()
            .map(|x| Ratio::<i128>::from_integer(x))
            .collect::<Vec<_>>();

        let series = comm![a, b];
        assert_eq!(series.coefficients, expected_coefficients);
    }
}
