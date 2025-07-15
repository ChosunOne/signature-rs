use std::ops::{Index, IndexMut};

use num_rational::Ratio;

use crate::{
    Commutator, Int,
    lyndon::{Generator, LyndonWord},
};

pub struct LieSeries<const N: usize, T: Generator<Letter = T>, U: Int + Send + Sync> {
    /// The Lyndon basis for the series
    basis: Vec<LyndonWord<N, T>>,
    /// The coefficients for each of the terms in the series
    coefficients: Vec<Ratio<U>>,
    /// The left factors of each word in the basis
    left_factors: Vec<usize>,
    /// The right factors of each word in the basis
    right_factors: Vec<usize>,
    /// The terms of the series
    terms: Vec<Ratio<U>>,
}

impl<const N: usize, T: Generator<Letter = T>, U: Int + Send + Sync> Index<usize>
    for LieSeries<N, T, U>
{
    type Output = Ratio<U>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.terms[index]
    }
}

impl<const N: usize, T: Generator<Letter = T>, U: Int + Send + Sync> IndexMut<usize>
    for LieSeries<N, T, U>
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.terms[index]
    }
}

impl<const N: usize, T: Generator<Letter = T>, U: Int + Send + Sync> LieSeries<N, T, U> {
    pub fn new(
        basis: Vec<LyndonWord<N, T>>,
        coefficients: Vec<Ratio<U>>,
        left_factors: Vec<usize>,
        right_factors: Vec<usize>,
        terms: Vec<Ratio<U>>,
    ) -> Self {
        Self {
            basis,
            coefficients,
            left_factors,
            right_factors,
            terms,
        }
    }
}

impl<const N: usize, T: Generator<Letter = T>, U: Int + Send + Sync> Commutator<&Self>
    for LieSeries<N, T, U>
{
    type Output = Self;

    fn commutator(&self, other: &Self) -> Self::Output {
        todo!()
    }
}
