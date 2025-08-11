use num_rational::Ratio;

use crate::{
    Int,
    lie_series::LieSeries,
    lyndon::{Generator, LyndonWord},
};
use std::{
    hash::Hash,
    ops::{Index, IndexMut},
};

pub struct LogSignature<const N: usize, T: Generator<Letter = T>, U: Int + Hash + Send + Sync> {
    series: LieSeries<N, T, U>,
    bch_series: LieSeries<N, T, U>,
}

impl<const N: usize, T: Generator<Letter = T>, U: Int + Hash + Send + Sync> Index<usize>
    for LogSignature<N, T, U>
{
    type Output = Ratio<U>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.series[index]
    }
}

impl<const N: usize, T: Generator<Letter = T>, U: Int + Hash + Send + Sync> Index<LyndonWord<N, T>>
    for LogSignature<N, T, U>
{
    type Output = Ratio<U>;

    fn index(&self, index: LyndonWord<N, T>) -> &Self::Output {
        &self.series[index]
    }
}

impl<const N: usize, T: Generator<Letter = T>, U: Int + Hash + Send + Sync> Index<&LyndonWord<N, T>>
    for LogSignature<N, T, U>
{
    type Output = Ratio<U>;

    fn index(&self, index: &LyndonWord<N, T>) -> &Self::Output {
        &self.series[index]
    }
}

impl<const N: usize, T: Generator<Letter = T>, U: Int + Hash + Send + Sync> IndexMut<usize>
    for LogSignature<N, T, U>
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.series[index]
    }
}

impl<const N: usize, T: Generator<Letter = T>, U: Int + Hash + Send + Sync>
    IndexMut<LyndonWord<N, T>> for LogSignature<N, T, U>
{
    fn index_mut(&mut self, index: LyndonWord<N, T>) -> &mut Self::Output {
        &mut self.series[index]
    }
}

impl<const N: usize, T: Generator<Letter = T>, U: Int + Hash + Send + Sync>
    IndexMut<&LyndonWord<N, T>> for LogSignature<N, T, U>
{
    fn index_mut(&mut self, index: &LyndonWord<N, T>) -> &mut Self::Output {
        &mut self.series[index]
    }
}

impl<const N: usize, T: Generator<Letter = T>, U: Int + Hash + Send + Sync> LogSignature<N, T, U> {
    pub fn new() -> Self {
        todo!()
    }

    pub fn concatenate(rhs: Self) -> Self {
        todo!()
    }
}
