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
    /// The coefficients for each of the terms in the series
    coefficients: Vec<Ratio<U>>,
    /// The left factors of each word in the basis
    left_factors: Vec<usize>,
    /// The right factors of each word in the basis
    right_factors: Vec<usize>,
    /// The terms of the series
    terms: Vec<Ratio<U>>,
}

impl<const N: usize, T: Generator<Letter = T>, U: Hash + Int + Send + Sync> Index<usize>
    for LieSeries<N, T, U>
{
    type Output = Ratio<U>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.terms[index]
    }
}

impl<const N: usize, T: Generator<Letter = T>, U: Hash + Int + Send + Sync> IndexMut<usize>
    for LieSeries<N, T, U>
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.terms[index]
    }
}

impl<const N: usize, T: Generator<Letter = T>, U: Hash + Int + Send + Sync> LieSeries<N, T, U> {
    pub fn new(
        basis: Vec<LyndonWord<N, T>>,
        coefficients: Vec<Ratio<U>>,
        left_factors: Vec<usize>,
        right_factors: Vec<usize>,
        terms: Vec<Ratio<U>>,
    ) -> Self {
        let mut commutator_basis = Vec::with_capacity(basis.len());
        for word in basis.iter() {
            commutator_basis.push(CommutatorTerm::from(word));
        }
        let mut commutator_basis_map = HashMap::new();

        for a in commutator_basis.iter() {
            for b in commutator_basis.iter() {
                if a == b {
                    continue;
                }
                let term = comm![a, b];
                let mut lyndon_term = term.clone();
                lyndon_term.lyndon_sort();
                commutator_basis_map.insert(term, lyndon_term);
            }
        }

        Self {
            basis,
            commutator_basis,
            commutator_basis_map,
            coefficients,
            left_factors,
            right_factors,
            terms,
        }
    }
}

impl<const N: usize, T: Generator<Letter = T> + Debug + Clone, U: Hash + Int + Send + Sync>
    Commutator<&Self> for LieSeries<N, T, U>
{
    type Output = Self;

    // For a series this is the following:
    // A: e_1 + 2e_2 + 3[e_1, e_2]
    // B: 4e_1 + 5e_2 + 6[e_1, e_2]
    //
    // [A, B] = [e_1, 4e_1] + [e_1, 5e_2] + [e_1, 6[e_1, e_2]]
    //          + [2e_2, 4e_1] + [2e_2, 5e_2] + [2e_2, 6[e_1, e_2]]
    //          + [3[e_1, e_2], 4e_1] + [3[e_1, e_2], 5e_2] + [3[e_1, e_2], 6[e_1, e_2]]
    //
    //        = [e_1, 5e_2] + [e_1, 6[e_1, e_2]]
    //          + [2e_2, 4e_1] + [2e_2, 6[e_1, e_2]]
    //          + [3[e_1, e_2], 4e_1] + [3[e_1, e_2], 5e_2]
    //
    //        = e_1*5e_2 - 5e_2*e_1 + e_1*6(e_1*e_2 - e_2*e_1) - 6(e_1*e_2 - e_2*e_1)*e_1
    //          + 2e_2*4e_1 - 4e_1*2e_2 + 2e_2 * 6(e_1*e_2 - e_2*e_1) - 6(e_1*e_2 - e_2*e_1)*2e_2
    //          + 3(e_1*e_2 - e_2*e_1)*4e_1 - 4e_1*3(e_1*e_2 - e_2*e_1)
    //          + 3(e_1*e_2 - e_2*e_1)*5e_2 - 5e_2*3(e_1*e_2 - e_2*e_1)
    //
    //        = 5e_1e_2 - 5e_2e_1 + 6e_1e_1e_2 - 6e_1e_2e_1 - 6e_1e_2e_1 + 6e_2e_1e_1
    //          + 8e_2e_1 - 8e_1e_2 + 12e_2e_1e_2 - 12e_2e_2e_1 - 12e_1e_2e_2 + 12e_2e_1e_2
    //          + 12e_1e_2e_1 - 12e_2e_1e_1 - 12e_1e_1e_2 + 12e_1e_2e_1
    //          + 15e_1e_2e_2 - 15e_2e_1e_2 - 15e_2e_1e_2 + 15e_2e_2e_1
    //
    //      Notational shortcut e_1e_2 = e_12
    //        = 5e_12 - 5e_21 + 6e_112 - 6e_121 - 6e_121 + 6e_211
    //          + 8e_21 - 8e_12 + 12e_212 - 12e_221 - 12e_122 + 12e_212
    //          + 12e_121 - 12e_211 - 12e_112 + 12e_121
    //          + 15e_122 - 15e_212 - 15e_212 + 15e_221
    //
    //        = 5e_12 - 8e_12
    //          + 8e_21 - 5e_21
    //          + 6e_112 - 12e_112
    //          + 12e_121 + 12e_121 - 6e_121 - 6e_121
    //          + 15e_122 - 12e_122
    //          + 6e_211 - 12e_211
    //          + 12e_212 + 12e_212 - 15e_212 - 15e_212
    //          + 15e_221 - 12e_221
    //
    //        = -3e_12 + 3e_21 - 6e_112 + 12e_121 + 3e_122 - 6e_211 - 6e_212 + 3e_221
    //
    //  Essentially, for two series A and B with bases:
    //  e_1 + e_2 + [e_1, e_2] + [e_1, [e_1, e_2]] + ...
    //
    //  we need to calculate all the possible combination terms:
    //  [e_1, e_1] + [e_1, e_2] + [e_1, [e_1, [e_1, e_2]]] + ...
    //  + [e_2, e_1] + [e_2, e_2] + [e_2, [e_1, e_2]] + ...
    //
    //  once we have that, we then need to store the corresponding lyndon basis element for each
    //  possible commutation.  Then we can use that lookup table to quickly aggregate the basis
    //  elements.

    fn commutator(&self, other: &Self) -> Self::Output {
        let mut terms = vec![U::default(); self.terms.len()];
        for i in 0..self.terms.len() {
            for j in 0..other.terms.len() {}
        }
        todo!()
    }
}
