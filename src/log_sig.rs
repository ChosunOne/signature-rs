use ndarray::{Array, Axis, Dimension, RemoveAxis};

use crate::{
    Arith, LieSeriesGenerator,
    bch_series_generator::BchSeriesGenerator,
    comm,
    commutator::{Commutator, CommutatorTerm},
    generators::Generator,
    lie_series::LieSeries,
    lyndon::{LyndonBasis, LyndonWord, Sort},
};
use ordered_float::NotNan;
use std::{
    collections::HashMap,
    fmt::Display,
    ops::{Index, IndexMut},
};

#[derive(Default, Debug)]
pub struct LogSignatureBuilder<T: Generator + Default + Send + Sync = u8> {
    /// The maximum degree of terms to include in the log signature
    max_degree: usize,
    lyndon_basis: LyndonBasis<T>,
}

impl<T: Generator + Default + Send + Sync> LogSignatureBuilder<T> {
    #[must_use]
    pub fn new() -> Self {
        Self {
            ..Default::default()
        }
    }

    #[must_use]
    pub fn with_max_degree(mut self, max_degree: usize) -> Self {
        self.max_degree = max_degree;
        self
    }

    #[must_use]
    pub fn with_num_dimensions(mut self, num_dimensions: usize) -> Self {
        self.lyndon_basis.alphabet_size = num_dimensions;
        self
    }

    #[must_use]
    pub fn max_degree(&self) -> usize {
        self.max_degree
    }

    #[must_use]
    pub fn num_dimensions(&self) -> usize {
        self.lyndon_basis.alphabet_size
    }

    #[must_use]
    pub fn build<U: Arith + Send + Sync>(&self) -> LogSignature<T, U> {
        let bch_basis = LyndonBasis::<u8>::new(2, Sort::Lexicographical);
        let bch_series = BchSeriesGenerator::new(bch_basis, self.max_degree).generate_lie_series();
        let basis = self.lyndon_basis.generate_basis(self.max_degree);
        let coefficients = vec![U::default(); basis.len()];
        let series = LieSeries::<T, U>::new(basis, coefficients);
        LogSignature::<T, U> { series, bch_series }
    }

    #[must_use]
    pub fn build_from_path<D: Dimension + RemoveAxis, U: Arith + Send + Sync>(
        &self,
        path: &Array<U, D>,
    ) -> LogSignature<T, U> {
        let mut log_sig = self.build();
        let mut log_sig_segment = log_sig.clone();

        for window in path.axis_windows(Axis(0), 2) {
            let start = window.index_axis(Axis(0), 0);
            let end = window.index_axis(Axis(0), 1);
            let displacement = (&end - &start).iter().cloned().collect::<Vec<_>>();
            log_sig_segment.series.coefficients[..self.lyndon_basis.alphabet_size]
                .clone_from_slice(&displacement[..self.lyndon_basis.alphabet_size]);
            log_sig.concatenate_assign(&log_sig_segment);
        }

        log_sig
    }
}

#[derive(Debug, Clone)]
pub struct LogSignature<T: Generator + Send + Sync = u8, U: Arith + Send + Sync = NotNan<f64>> {
    pub series: LieSeries<T, U>,
    pub bch_series: LieSeries<u8, U>,
}

impl<T: Generator + Send + Sync, U: Arith + Display + Send + Sync> Index<usize>
    for LogSignature<T, U>
{
    type Output = U;

    fn index(&self, index: usize) -> &Self::Output {
        &self.series[index]
    }
}

impl<T: Generator + Send + Sync, U: Arith + Display + Send + Sync> Index<LyndonWord<T>>
    for LogSignature<T, U>
{
    type Output = U;

    fn index(&self, index: LyndonWord<T>) -> &Self::Output {
        &self.series[index]
    }
}

impl<T: Generator + Send + Sync, U: Arith + Display + Send + Sync> Index<&LyndonWord<T>>
    for LogSignature<T, U>
{
    type Output = U;

    fn index(&self, index: &LyndonWord<T>) -> &Self::Output {
        &self.series[index]
    }
}

impl<T: Generator + Send + Sync, U: Arith + Display + Send + Sync> IndexMut<usize>
    for LogSignature<T, U>
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.series[index]
    }
}

impl<T: Generator + Send + Sync, U: Arith + Display + Send + Sync> IndexMut<LyndonWord<T>>
    for LogSignature<T, U>
{
    fn index_mut(&mut self, index: LyndonWord<T>) -> &mut Self::Output {
        &mut self.series[index]
    }
}

impl<T: Generator + Send + Sync, U: Arith + Display + Send + Sync> IndexMut<&LyndonWord<T>>
    for LogSignature<T, U>
{
    fn index_mut(&mut self, index: &LyndonWord<T>) -> &mut Self::Output {
        &mut self.series[index]
    }
}

impl<T: Generator + Send + Sync, U: Arith + Send + Sync> LogSignature<T, U> {
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

    pub fn concatenate_assign(&mut self, rhs: &Self) {
        let mut computed_commutations = HashMap::new();
        let original_series = self.series.clone();

        for (i, term) in self.bch_series.commutator_basis.iter().enumerate().skip(1) {
            self.series += evaluate_commutator_term(
                term,
                &[&original_series, &rhs.series],
                &mut computed_commutations,
            ) * self.bch_series[i].clone();
        }
    }
}

fn evaluate_commutator_term<T: Generator, U: Arith + Send + Sync>(
    term: &CommutatorTerm<U, u8>,
    series: &[&LieSeries<T, U>],
    computed_commutations: &mut HashMap<CommutatorTerm<U, u8>, LieSeries<T, U>>,
) -> LieSeries<T, U> {
    match term {
        &CommutatorTerm::Atom { atom: a, .. } => series[a as usize].clone(),
        t @ CommutatorTerm::Expression { left, right, .. } => {
            if computed_commutations.contains_key(t) {
                return computed_commutations[t].clone();
            }
            let a = evaluate_commutator_term(left, series, computed_commutations);
            let b = evaluate_commutator_term(right, series, computed_commutations);
            let result = comm![a, b];
            computed_commutations.insert(t.clone(), result.clone());
            result
        }
    }
}

#[cfg(test)]
mod test {
    use ndarray::{Array2, array};
    use num_rational::Ratio;
    use num_traits::ToPrimitive;
    use rstest::rstest;

    use crate::generators::ENotation;

    use super::*;

    #[rstest]
    #[case(
        3,
        3,
        array![
            [0.0, 0., 0.],
            [1.0, 2.0, 3.0],
        ],
        vec![
            1.000000,
            2.000000,
            3.000000,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
        ])]
    #[case(
        3,
        3,
        array![
            [0.0, 0., 0.],
            [1.0, 2.0, 3.0],
            [6.0, 5.0, 4.0],
        ],
        vec![
            6.,
            5.,
            4.,
            -3.5,
            -7.,
            -3.5,
            2.333333,
            4.666667,
            -0.583333,
            3.5,
            0.,
            2.333333,
            0.583333,
            1.166667,
        ])]
    #[case(
        3,
        3,
        array![
            [0.0, 0., 0.],
            [1.0, 2.0, 3.0],
            [6.0, 5.0, 4.0],
            [7.0, 8.0, 9.0],
            [12.0, 11.0, 10.0]
        ],
        vec![
            12.000000,
            11.000000,
            10.000000,
            -6.500000,
            -13.000000,
            -6.500000,
            -1.166667,
            -2.333333,
            4.416667,
            6.500000,
            16.500000,
            15.333333,
            -4.416667,
            7.666667,
        ])]
    #[case(
        3,
        3,
        array![
            [0.0, 0., 0.],
            [-0.077, 0.042, -0.067],
            [-0.154, 0.675, 0.006],
            [0.916, 1.177, -0.139],
            [1.095, 0.823, -0.261]
        ],
        vec![
           1.095000,
           0.823000,
          -0.261000,
          -0.690006,
          -0.040871,
          -0.124105,
           0.098690,
          -0.004304,
           0.146613,
           0.024960,
           0.044713,
          -0.000215,
          -0.038903,
           0.001568,
        ])]
    #[case(
        3,
        4,
        array![
            [0., 0., 0.],
            [1.0, 2.0, 3.0],
            [6.0, 5.0, 4.0],
        ],
        vec![
            6.000000,
            5.000000,
            4.000000,
            -3.500000,
            -7.000000,
            -3.500000,
            2.333333,
            4.666667,
            -0.583333,
            3.500000,
            0.000000,
            2.333333,
            0.583333,
            1.166667,
            1.458333,
            2.916667,
            -3.791667,
            -3.208333,
            -12.250000,
            -9.333333,
            -1.458333,
            1.750000,
            3.500000,
            2.916667,
            -3.791667,
            6.708333,
            2.625000,
            7.291667,
            1.750000,
            1.750000,
            -3.208333,
            0.875000,
        ]
    )]
    #[case(
        3,
        4,
        array![
            [   0.000,    0.000,    0.000],
            [   1.000,    2.000,    3.000],
            [   6.000,    5.000,    4.000],
            [   7.000,    8.000,    9.000],
            [  12.000,   11.000,   10.000]
        ],
        vec![
            12.000000,
            11.000000,
            10.000000,
            -6.500000,
            -13.000000,
            -6.500000,
            -1.166667,
            -2.333333,
            4.416667,
            6.500000,
            16.500000,
            15.333333,
            -4.416667,
            7.666667,
            -0.041667,
            -0.083333,
            1.208333,
            2.291667,
            4.750000,
            4.666667,
            0.041667,
            -2.250000,
            -4.500000,
            -11.083333,
            -4.291667,
            -12.291667,
            -19.875000,
            -22.208333,
            -13.250000,
            -2.250000,
            7.791667,
            -6.625000,
        ]
    )]
    #[case(
        4,
        4,
        array![
            [   0.000,    0.000,    0.000,    0.000],
            [   1.000,    2.000,    3.000,    4.000],
            [   6.000,    5.000,    4.000,    3.000],
            [   7.000,    8.000,    9.000,    8.000],
            [  12.000,   11.000,   10.000,    9.000],
        ],
        vec![
    12.000000,
      11.000000,
      10.000000,
       9.000000,
      -6.500000,
     -13.000000,
     -13.500000,
      -6.500000,
      -7.000000,
      -0.500000,
      -1.166667,
      -2.333333,
      10.500000,
       4.416667,
       6.500000,
      18.583333,
      16.500000,
      15.333333,
      26.666667,
       5.166667,
      20.833333,
       7.750000,
      -4.416667,
      6.166667,
       7.666667,
      17.500000,
       6.250000,
       0.833333,
       8.333333,
      -6.083333,
      -0.041667,
      -0.083333,
     -14.625000,
       1.208333,
       2.291667,
     -19.125000,
       4.750000,
       4.666667,
     -23.625000,
      21.083333,
      12.916667,
       2.375000,
       0.041667,
       8.083333,
      -2.250000,
      -4.500000,
     -18.750000,
     -11.083333,
      -4.291667,
     -24.500000,
       4.583333,
       8.416667,
       5.750000,
       1.541667,
     -12.291667,
     -19.875000,
      -9.958333,
     -22.208333,
     -13.250000,
     -25.375000,
       1.083333,
     -21.458333,
       9.125000,
     -14.833333,
     -16.666667,
      -9.500000,
     -38.250000,
     -31.791667,
     -19.000000,
     -12.416667,
     -22.458333,
      -4.500000,
      -2.250000,
     -11.500000,
       7.791667,
      -8.166667,
      22.666667,
      10.166667,
       4.250000,
      -6.625000,
     -15.250000,
      -8.166667,
       7.916667,
     -18.458333,
      -9.500000,
     -14.583333,
      -3.333333,
      -5.125000,
       6.708333,
      -2.166667,
        ]
    )]
    fn test_log_sig_builder_from_path(
        #[case] num_dimensions: usize,
        #[case] max_degree: usize,
        #[case] path: Array2<f64>,
        #[case] expected_coefficients: Vec<f64>,
    ) {
        let builder = LogSignatureBuilder::<ENotation>::new()
            .with_max_degree(max_degree)
            .with_num_dimensions(num_dimensions);
        let path = path.mapv(|v| NotNan::new(v).expect("value not to be a number"));
        let log_sig = builder.build_from_path(&path);
        for (i, c) in log_sig.series.commutator_basis.iter().enumerate() {
            println!("{i}: {} \t {c}", log_sig.series.basis[i]);
        }
        for (i, c) in log_sig.series.coefficients.iter().enumerate() {
            println!("[{i}]: {c}");
        }
        assert_eq!(
            log_sig.series.coefficients.len(),
            expected_coefficients.len()
        );
        for (i, &c) in expected_coefficients.iter().enumerate() {
            assert!(
                (c - log_sig.series.coefficients[i].to_f64().unwrap()).abs() < 0.0001,
                "{i}: {} != {c}",
                log_sig.series.coefficients[i].to_f64().unwrap()
            );
        }
    }

    #[test]
    fn test_log_sig_concat() {
        let builder = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(2)
            .with_max_degree(3);
        let mut a = builder.build();
        let mut b = builder.build();
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
