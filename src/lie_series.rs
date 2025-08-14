use std::collections::HashMap;
use std::fmt::Debug;
use std::ops::{Add, AddAssign, Index, IndexMut, Mul, Sub, SubAssign};

use crate::Arith;
use crate::{
    comm,
    commutator::{Commutator, CommutatorTerm},
    lyndon::{Generator, LyndonWord},
};

#[derive(Debug, Default, Clone)]
pub struct LieSeries<T: Generator, U: Arith + Send + Sync> {
    /// The Lyndon basis for the series
    basis: Vec<LyndonWord<T>>,
    /// The commutator basis for the series
    pub commutator_basis: Vec<CommutatorTerm<U, T>>,
    /// A map for converting arbitrary commutator terms to basis elements
    commutator_basis_map: HashMap<CommutatorTerm<U, T>, CommutatorTerm<U, T>>,
    /// A map for locating a given term's index in the basis
    commutator_basis_index_map: HashMap<CommutatorTerm<U, T>, usize>,
    /// The coefficients for each of the terms in the series
    pub coefficients: Vec<U>,
}

impl<T: Generator, U: Arith + Send + Sync> Index<usize> for LieSeries<T, U> {
    type Output = U;

    fn index(&self, index: usize) -> &Self::Output {
        &self.coefficients[index]
    }
}

impl<T: Generator, U: Arith + Send + Sync> Index<LyndonWord<T>> for LieSeries<T, U> {
    type Output = U;

    fn index(&self, index: LyndonWord<T>) -> &Self::Output {
        let term = CommutatorTerm::<U, T>::from(&index);
        let i = self.commutator_basis_index_map[&term];
        &self.coefficients[i]
    }
}

impl<T: Generator, U: Arith + Send + Sync> Index<&LyndonWord<T>> for LieSeries<T, U> {
    type Output = U;

    fn index(&self, index: &LyndonWord<T>) -> &Self::Output {
        let term = CommutatorTerm::<U, T>::from(index);
        let i = self.commutator_basis_index_map[&term];
        &self.coefficients[i]
    }
}

impl<T: Generator, U: Arith + Send + Sync> IndexMut<usize> for LieSeries<T, U> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.coefficients[index]
    }
}

impl<T: Generator, U: Arith + Send + Sync> IndexMut<LyndonWord<T>> for LieSeries<T, U> {
    fn index_mut(&mut self, index: LyndonWord<T>) -> &mut Self::Output {
        let term = CommutatorTerm::<U, T>::from(&index);
        let i = self.commutator_basis_index_map[&term];
        &mut self.coefficients[i]
    }
}

impl<T: Generator, U: Arith + Send + Sync> IndexMut<&LyndonWord<T>> for LieSeries<T, U> {
    fn index_mut(&mut self, index: &LyndonWord<T>) -> &mut Self::Output {
        let term = CommutatorTerm::<U, T>::from(index);
        let i = self.commutator_basis_index_map[&term];
        &mut self.coefficients[i]
    }
}

impl<T: Generator, U: Arith + Clone + Send + Sync> Add for &LieSeries<T, U> {
    type Output = LieSeries<T, U>;

    fn add(self, rhs: Self) -> Self::Output {
        let mut coefficients = self.coefficients.clone();
        for i in 0..self.coefficients.len() {
            coefficients[i] = self.coefficients[i].clone() + rhs.coefficients[i].clone();
        }
        LieSeries::<T, U> {
            basis: self.basis.clone(),
            commutator_basis: self.commutator_basis.clone(),
            commutator_basis_map: self.commutator_basis_map.clone(),
            commutator_basis_index_map: self.commutator_basis_index_map.clone(),
            coefficients,
        }
    }
}

impl<T: Generator, U: Arith + Send + Sync> Add for LieSeries<T, U> {
    type Output = LieSeries<T, U>;

    fn add(mut self, rhs: Self) -> Self::Output {
        for i in 0..self.coefficients.len() {
            self.coefficients[i] += rhs.coefficients[i].clone();
        }
        self
    }
}

impl<T: Generator, U: Arith + Send + Sync> AddAssign for LieSeries<T, U> {
    fn add_assign(&mut self, rhs: Self) {
        for i in 0..self.coefficients.len() {
            self.coefficients[i] += rhs.coefficients[i].clone();
        }
    }
}

impl<T: Generator, U: Arith + Send + Sync> AddAssign<&Self> for LieSeries<T, U> {
    fn add_assign(&mut self, rhs: &Self) {
        for i in 0..self.coefficients.len() {
            self.coefficients[i] += rhs.coefficients[i].clone();
        }
    }
}

impl<T: Generator, U: Arith + Send + Sync> Sub for &LieSeries<T, U> {
    type Output = LieSeries<T, U>;

    fn sub(self, rhs: Self) -> Self::Output {
        let mut coefficients = self.coefficients.clone();
        for i in 0..self.coefficients.len() {
            coefficients[i] = self.coefficients[i].clone() - rhs.coefficients[i].clone();
        }
        LieSeries::<T, U> {
            basis: self.basis.clone(),
            commutator_basis: self.commutator_basis.clone(),
            commutator_basis_map: self.commutator_basis_map.clone(),
            commutator_basis_index_map: self.commutator_basis_index_map.clone(),
            coefficients,
        }
    }
}

impl<T: Generator, U: Arith + Send + Sync> Sub for LieSeries<T, U> {
    type Output = LieSeries<T, U>;

    fn sub(mut self, rhs: Self) -> Self::Output {
        for i in 0..self.coefficients.len() {
            self.coefficients[i] -= rhs.coefficients[i].clone();
        }
        self
    }
}

impl<T: Generator, U: Arith + Send + Sync> SubAssign for LieSeries<T, U> {
    fn sub_assign(&mut self, rhs: Self) {
        for i in 0..self.coefficients.len() {
            self.coefficients[i] -= rhs.coefficients[i].clone();
        }
    }
}

impl<T: Generator, U: Arith + Send + Sync> SubAssign<&Self> for LieSeries<T, U> {
    fn sub_assign(&mut self, rhs: &Self) {
        for i in 0..self.coefficients.len() {
            self.coefficients[i] -= rhs.coefficients[i].clone();
        }
    }
}

impl<T: Generator, U: Arith + Send + Sync> Mul<U> for LieSeries<T, U> {
    type Output = Self;

    fn mul(mut self, rhs: U) -> Self::Output {
        for i in 0..self.coefficients.len() {
            self.coefficients[i] *= rhs.clone();
        }
        self
    }
}

impl<T: Generator, U: Arith + Send + Sync> Mul<U> for &LieSeries<T, U> {
    type Output = LieSeries<T, U>;

    fn mul(self, rhs: U) -> Self::Output {
        let mut result = self.clone();
        for i in 0..self.coefficients.len() {
            result.coefficients[i] *= rhs.clone();
        }
        result
    }
}

impl<T: Generator, U: Arith + Send + Sync> LieSeries<T, U> {
    #[must_use]
    pub fn new(basis: Vec<LyndonWord<T>>, coefficients: Vec<U>) -> Self {
        let mut commutator_basis = Vec::<CommutatorTerm<U, T>>::with_capacity(basis.len());
        for word in &basis {
            commutator_basis.push(CommutatorTerm::from(word));
        }
        let mut commutator_basis_map = HashMap::new();

        let mut commutator_basis_index_map = HashMap::new();
        for (i, a) in commutator_basis.iter().enumerate() {
            commutator_basis_index_map.insert(a.clone(), i);
            if let CommutatorTerm::Expression(e) = a {
                let mut e = e.clone();
                e.coefficient = -e.coefficient.clone();
                commutator_basis_index_map.insert(CommutatorTerm::Expression(e), i);
            }
        }

        for a in &commutator_basis {
            for b in &commutator_basis {
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

impl<T: Generator + Debug + Clone, U: Arith + Send + Sync> Commutator<&Self> for LieSeries<T, U> {
    type Output = Self;

    /// Calculates the lie bracket `[A, B]` for a lie series for terms within the commutator basis.
    fn commutator(&self, other: &Self) -> Self::Output {
        let mut coefficients = vec![U::default(); self.coefficients.len()];
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
        use num_rational::Ratio;

        use crate::lyndon::{LyndonBasis, Sort};

        let basis = LyndonBasis::<u8>::new(2, Sort::Lexicographical).generate_basis(basis_depth);
        let a_coefficients = a_coefficients
            .into_iter()
            .map(Ratio::<i128>::from_integer)
            .collect::<Vec<_>>();
        let a = LieSeries::new(basis.clone(), a_coefficients);
        let b_coefficients = b_coefficients
            .into_iter()
            .map(Ratio::<i128>::from_integer)
            .collect::<Vec<_>>();
        let b = LieSeries::new(basis.clone(), b_coefficients);
        let expected_coefficients = expected_coefficients
            .into_iter()
            .map(Ratio::<i128>::from_integer)
            .collect::<Vec<_>>();

        let series = comm![a, b];
        assert_eq!(series.coefficients, expected_coefficients);
    }
}
