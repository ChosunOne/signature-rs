use lie_rs::LieSeriesGenerator;
use lie_rs::{LieSeries, bch_series_generator::BchSeriesGenerator};
use lyndon_rs::lyndon::LyndonWord;
use lyndon_rs::{
    generators::Generator,
    lyndon::{LyndonBasis, Sort},
};
use ndarray::{ArrayView, Axis, Dimension, RemoveAxis};

use num_traits::{FromPrimitive, One, Zero};
use std::ops::{Mul, Sub};
use std::{
    fmt::Debug,
    hash::Hash,
    ops::{AddAssign, Div, Index, IndexMut, MulAssign, Neg, SubAssign},
};

use crate::commutator_dag::{CommutatorDag, TermSource};

/// Builder for constructing log signatures from path data.
///
/// The log signature is a mathematical transform that captures the geometry
/// of a path by expressing it as a series in the free Lie algebra. This builder
/// allows configuring the parameters and constructing log signatures from various inputs.
pub struct LogSignatureBuilder<T> {
    /// The maximum degree of terms to include in the log signature computation.
    pub max_degree: usize,
    /// The Lyndon basis configuration for the underlying algebra.
    pub lyndon_basis: LyndonBasis<T>,
}

impl<T: Debug> Debug for LogSignatureBuilder<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("LogSignatureBuilder")
            .field("max_degree", &self.max_degree)
            .field("lyndon_basis", &self.lyndon_basis)
            .finish()
    }
}

impl<T> Copy for LogSignatureBuilder<T> {}

impl<T> Clone for LogSignatureBuilder<T> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<T> Default for LogSignatureBuilder<T> {
    fn default() -> Self {
        Self {
            max_degree: usize::default(),
            lyndon_basis: LyndonBasis::default(),
        }
    }
}

impl<T> LogSignatureBuilder<T> {
    /// Creates a new log signature builder with default settings.
    #[must_use]
    pub fn new() -> Self {
        Self {
            ..Default::default()
        }
    }

    /// Sets the maximum degree of terms to include in the log signature.
    ///
    /// Higher degrees capture more complex geometric features but increase
    /// computational complexity exponentially.
    #[must_use]
    pub fn with_max_degree(mut self, max_degree: usize) -> Self {
        self.max_degree = max_degree;
        self
    }

    /// Sets the number of dimensions for the path data.
    ///
    /// This determines the size of the generator alphabet and should match
    /// the dimensionality of the input path data.
    #[must_use]
    pub fn with_num_dimensions(mut self, num_dimensions: usize) -> Self {
        self.lyndon_basis.alphabet_size = num_dimensions;
        self
    }

    /// Returns the maximum degree setting.
    #[must_use]
    pub fn max_degree(&self) -> usize {
        self.max_degree
    }

    /// Returns the number of dimensions setting.
    #[must_use]
    pub fn num_dimensions(&self) -> usize {
        self.lyndon_basis.alphabet_size
    }
}

impl<T: Debug + Clone + Eq + Hash + Ord + Generator + Send + Sync> LogSignatureBuilder<T> {
    /// Builds an empty log signature with the configured parameters.
    ///
    /// The resulting log signature has the proper basis structure but with
    /// all coefficients set to zero.
    #[must_use]
    pub fn build<
        U: Clone
            + Default
            + AddAssign
            + Div<Output = U>
            + FromPrimitive
            + Neg<Output = U>
            + One
            + Ord
            + Hash
            + Send
            + Sync
            + Zero
            + SubAssign
            + MulAssign
            + Send
            + Sync,
    >(
        &self,
    ) -> LogSignature<T, U> {
        let bch_basis = LyndonBasis::<u8>::new(2, Sort::Lexicographical);
        let bch_series = BchSeriesGenerator::new(bch_basis, self.max_degree).generate_lie_series();
        let dag = CommutatorDag::from_bch_series(&bch_series);
        let basis = self.lyndon_basis.generate_basis(self.max_degree);
        let coefficients = vec![U::default(); basis.len()];
        let series = LieSeries::<T, U>::new(basis, coefficients);
        LogSignature::<T, U> {
            series,
            bch_series,
            dag,
        }
    }

    /// Computes the log signature of a path from multidimensional array data.
    ///
    /// The path should be provided as a 2D array where each row represents a point
    /// and each column represents a coordinate dimension. The log signature is
    /// computed incrementally over consecutive path segments.
    #[must_use]
    pub fn build_from_path<
        D: Dimension + RemoveAxis,
        U: Clone
            + Default
            + AddAssign
            + Div<Output = U>
            + FromPrimitive
            + Neg<Output = U>
            + One
            + Ord
            + Hash
            + Send
            + Sync
            + Zero
            + SubAssign
            + MulAssign
            + Sub<Output = U>,
    >(
        &self,
        path: &ArrayView<U, D>,
    ) -> LogSignature<T, U> {
        let mut log_sig = self.build();

        for window in path.axis_windows(Axis(0), 2) {
            let start = window.index_axis(Axis(0), 0);
            let end = window.index_axis(Axis(0), 1);
            let displacement = (&end - &start).iter().cloned().collect::<Vec<_>>();
            log_sig.concatenate_assign_coefficients(&displacement);
        }

        log_sig
    }
}

/// Represents a computed log signature of a path.
///
/// A log signature captures the essential geometric and algebraic features of a path
/// through its representation as a series in the free Lie algebra. This structure
/// contains both the computed coefficients and the Baker-Campbell-Hausdorff series
/// needed for concatenation operations.
pub struct LogSignature<T, U> {
    /// The main series containing the log signature coefficients.
    pub series: LieSeries<T, U>,
    /// The BCH series used for concatenating log signatures.
    pub bch_series: LieSeries<u8, U>,
    /// Compiled evaluation plan for [`LogSignature::concatenate_assign`]
    dag: CommutatorDag<U>,
}

impl<T: Debug, U: Debug> Debug for LogSignature<T, U> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("LogSignature")
            .field("series", &self.series)
            .field("bch_series", &self.bch_series)
            .field("dag", &self.dag)
            .finish()
    }
}

impl<T: Clone, U: Clone> Clone for LogSignature<T, U> {
    fn clone(&self) -> Self {
        Self {
            series: self.series.clone(),
            bch_series: self.bch_series.clone(),
            dag: self.dag.clone_shallow(),
        }
    }
}

impl<T, U> Index<usize> for LogSignature<T, U> {
    type Output = U;

    fn index(&self, index: usize) -> &Self::Output {
        &self.series[index]
    }
}

impl<
    T: Clone + Ord + Generator + Hash,
    U: Clone + One + Zero + Eq + MulAssign + Neg<Output = U> + Hash,
> Index<LyndonWord<T>> for LogSignature<T, U>
{
    type Output = U;

    fn index(&self, index: LyndonWord<T>) -> &Self::Output {
        &self.series[index]
    }
}

impl<
    T: Clone + Ord + Generator + Hash,
    U: Clone + One + Zero + Eq + MulAssign + Neg<Output = U> + Hash,
> Index<&LyndonWord<T>> for LogSignature<T, U>
{
    type Output = U;

    fn index(&self, index: &LyndonWord<T>) -> &Self::Output {
        &self.series[index]
    }
}

impl<T, U> IndexMut<usize> for LogSignature<T, U> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.series[index]
    }
}

impl<
    T: Clone + Ord + Generator + Hash,
    U: Clone + One + Zero + Eq + MulAssign + Neg<Output = U> + Hash,
> IndexMut<LyndonWord<T>> for LogSignature<T, U>
{
    fn index_mut(&mut self, index: LyndonWord<T>) -> &mut Self::Output {
        &mut self.series[index]
    }
}

impl<
    T: Clone + Ord + Generator + Hash,
    U: Clone + One + Zero + Eq + MulAssign + Neg<Output = U> + Hash,
> IndexMut<&LyndonWord<T>> for LogSignature<T, U>
{
    fn index_mut(&mut self, index: &LyndonWord<T>) -> &mut Self::Output {
        &mut self.series[index]
    }
}

impl<
    T: Clone + Ord + Generator + Hash,
    U: Clone
        + Mul<Output = U>
        + MulAssign
        + AddAssign
        + Hash
        + Eq
        + Default
        + One
        + Zero
        + Neg<Output = U>,
> LogSignature<T, U>
{
    /// Concatenates two log signatures using the Baker-Campbell-Hausdorff formula.
    ///
    /// This operation computes the log signature of the path formed by concatenating
    /// the paths represented by `self` and `rhs`. The result captures the geometry
    /// of the combined path.
    #[must_use]
    pub fn concatenate(&self, rhs: &Self) -> Self
    where
        U: Send + Sync,
    {
        let mut concatenated_log_sig = self.clone();

        concatenated_log_sig.concatenate_assign(rhs);

        concatenated_log_sig
    }

    /// In-place concatenation of log signatures.
    ///
    /// This is the mutable version of `concatenate` that modifies `self` instead
    /// of creating a new log signature. More memory-efficient for chaining operations.
    pub fn concatenate_assign(&mut self, rhs: &Self)
    where
        U: Send + Sync,
    {
        self.concatenate_assign_coefficients(&rhs.series.coefficients);
    }

    pub(crate) fn concatenate_assign_coefficients(&mut self, rhs_coefficients: &[U])
    where
        U: Send + Sync,
    {
        let original_coefficients = self.series.coefficients.clone();

        let a_nonzero = self
            .series
            .nonzero_coefficient_indices(&original_coefficients);
        let b_nonzero = self.series.nonzero_coefficient_indices(rhs_coefficients);

        self.dag.evaluate(
            &self.series,
            &original_coefficients,
            &a_nonzero,
            rhs_coefficients,
            &b_nonzero,
        );

        for (source, weight) in self.dag.terms() {
            let ct: &[U] = match source {
                TermSource::Node(node) => self.dag.buffer(node),
                TermSource::Displacement => rhs_coefficients,
            };
            for (self_coeff, comm_coeff) in self.series.coefficients.iter_mut().zip(ct.iter()) {
                *self_coeff += comm_coeff.clone() * weight.clone();
            }
        }
    }
}

#[cfg(test)]
mod test {
    use ndarray::{Array2, array};
    use num_rational::Ratio;
    use num_traits::ToPrimitive;
    use rstest::rstest;

    use super::*;

    #[test]
    fn dag_node_lists_cover_buffer_support() {
        use ordered_float::NotNan;

        let (d, m) = (3usize, 5usize);
        let builder: LogSignatureBuilder<u8> = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(d)
            .with_max_degree(m);
        let mut log_sig: LogSignature<u8, NotNan<f64>> = builder.build::<NotNan<f64>>();

        // Dense accumulator, then several letter-displacement folds: the   Not Committed Yet
        // first fold goes through the collecting rebuild, later folds reuse
        // the fixed-point lists.
        let basis: Vec<LyndonWord<u8>> =
            lyndon_rs::lyndon::LyndonBasis::<u8>::new(d, lyndon_rs::lyndon::Sort::Lexicographical)
                .generate_basis(m);
        let mut seed: u64 = 0xabcd_u64;
        let mut lcg = |seed: &mut u64| {
            *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((*seed >> 33) as i64) as f64
        };
        let acc: Vec<NotNan<f64>> = (0..basis.len())
            .map(|_| NotNan::new(lcg(&mut seed)).unwrap())
            .collect();
        log_sig.series.coefficients.clone_from(&acc);

        for _ in 0..4 {
            let disp: Vec<NotNan<f64>> = (0..basis.len())
                .map(|k| NotNan::new(if k < d { lcg(&mut seed) } else { 0.0 }).unwrap())
                .collect();
            let mut disp_sig: LogSignature<u8, NotNan<f64>> = builder.build::<NotNan<f64>>();
            disp_sig.series.coefficients.clone_from(&disp);
            log_sig.concatenate_assign(&disp_sig);

            let dag: &CommutatorDag<NotNan<f64>> = &log_sig.dag;
            let cutoff: usize = log_sig
                .series
                .commutator_basis
                .iter()
                .take_while(|w| w.degree() != m)
                .count();
            for k in 2..dag.structure.nodes.len() {
                let buffer: &Vec<NotNan<f64>> = &dag.buffers[k - 2];
                let list: &Vec<usize> = &dag.nonzeros[k - 2];
                let mut sorted: Vec<usize> = list.clone();
                sorted.sort_unstable();
                sorted.dedup();
                assert_eq!(
                    sorted.len(),
                    list.len(),
                    "node {k}: list contains duplicate indices"
                );
                for (i, c) in buffer.iter().enumerate().take(cutoff) {
                    assert!(
                        c.is_zero() || list.contains(&i),
                        "node {k}: index {i} is non-zero but not listed"
                    );
                }
            }
        }
    }

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
        use lyndon_rs::generators::ENotation;
        use ndarray::s;
        use ordered_float::NotNan;

        let builder = LogSignatureBuilder::<ENotation>::new()
            .with_max_degree(max_degree)
            .with_num_dimensions(num_dimensions);
        let path = path.mapv(|v| NotNan::new(v).expect("value to be a number"));
        let log_sig = builder.build_from_path(&path.slice(s![.., ..]));
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

    #[rstest]
    #[case(5, 5, array![
        [1., 2., 3., 4., 5.],
        [6., 7., 8., 9., 10.],
        [11., 12., 13., 14., 15.],
        [16., 17., 18., 19., 20.],
    ])]
    fn test_log_sig_concat_from_path(
        #[case] num_dimensions: usize,
        #[case] max_degree: usize,
        #[case] path: Array2<f64>,
    ) {
        use ndarray::s;
        use ordered_float::NotNan;

        let builder = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(num_dimensions)
            .with_max_degree(max_degree);
        let path = path.mapv(|v| NotNan::new(v).expect("value to be a number"));

        dbg!(&path.slice(s![0..2, ..]));

        let mut concatenated_log_sig = builder.build_from_path(&path.slice(s![0..=1, ..]));
        let log_sig_2 = builder.build_from_path(&path.slice(s![1.., ..]));

        concatenated_log_sig.concatenate_assign(&log_sig_2);

        let log_sig = builder.build_from_path(&path.slice(s![.., ..]));

        for (concat_c, full_path_c) in concatenated_log_sig
            .series
            .coefficients
            .iter()
            .zip(log_sig.series.coefficients.iter())
        {
            assert_eq!(concat_c, full_path_c);
        }
    }

    /// Differential test: after the first fold (which builds the node lists
    /// through the collecting pass), the steady level-batch evaluation must
    /// produce bit-identical node buffers to a fresh DAG's collecting
    /// rebuild — the oracle — for every subsequent fold.
    #[test]
    fn steady_batch_matches_rebuild_oracle() {
        use ordered_float::NotNan;

        for (d, m, folds) in [(3usize, 3usize, 6), (3, 4, 6), (2, 8, 6), (4, 5, 4)] {
            let builder = LogSignatureBuilder::<u8>::new()
                .with_num_dimensions(d)
                .with_max_degree(m);
            let mut log_sig = builder.build::<NotNan<f64>>();

            let basis = lyndon_rs::lyndon::LyndonBasis::<u8>::new(
                d,
                lyndon_rs::lyndon::Sort::Lexicographical,
            )
            .generate_basis(m);
            let mut seed = 0x5eed_u64;
            let mut lcg = |seed: &mut u64| {
                *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                ((*seed >> 33) as i64) as f64
            };
            let acc: Vec<NotNan<f64>> = (0..basis.len())
                .map(|_| NotNan::new(lcg(&mut seed)).unwrap())
                .collect();
            log_sig.series.coefficients.clone_from(&acc);

            for fold in 0..folds {
                // Alternate letter-only displacements (the production fold's
                // regime) with dense ones (worst case for the gating).
                let dense = fold % 2 == 1;
                let disp: Vec<NotNan<f64>> = (0..basis.len())
                    .map(|k| {
                        NotNan::new(if dense || k < d { lcg(&mut seed) } else { 0.0 }).unwrap()
                    })
                    .collect();
                let mut disp_sig = builder.build::<NotNan<f64>>();
                disp_sig.series.coefficients.clone_from(&disp);

                let original = log_sig.series.coefficients.clone();
                let a_nonzero = log_sig.series.nonzero_coefficient_indices(&original);
                let b_nonzero = log_sig.series.nonzero_coefficient_indices(&disp);
                log_sig.dag.evaluate(
                    &log_sig.series,
                    &original,
                    &a_nonzero,
                    &disp,
                    &b_nonzero,
                );

                // Oracle: a fresh DAG has no lists, so its evaluate takes the
                // collecting-rebuild path with the identical inputs.
                let mut oracle = CommutatorDag::<NotNan<f64>>::from_bch_series(
                    &log_sig.bch_series,
                );
                oracle.evaluate(&log_sig.series, &original, &a_nonzero, &disp, &b_nonzero);
                for (k, (buf, reff)) in
                    log_sig.dag.buffers.iter().zip(&oracle.buffers).enumerate()
                {
                    assert_eq!(
                        buf.len(),
                        reff.len(),
                        "d={d} m={m} fold={fold}: node {} buffer lengths differ",
                        k + 2
                    );
                    for (i, (x, y)) in buf.iter().zip(reff.iter()).enumerate() {
                        assert!(
                            x == y,
                            "d={d} m={m} fold={fold}: node {} index {} differs: {x:?} vs {y:?}",
                            k + 2,
                            i
                        );
                    }
                }

                // Complete the fold so the next one sees the evolved
                // accumulator (and the steady-batch lists).
                for (source, weight) in log_sig.dag.terms() {
                    let ct: &[NotNan<f64>] = match source {
                        TermSource::Node(node) => log_sig.dag.buffer(node),
                        TermSource::Displacement => &disp,
                    };
                    for (c, comm) in log_sig.series.coefficients.iter_mut().zip(ct) {
                        *c += comm.clone() * weight.clone();
                    }
                }
            }
        }
    }

    /// Fold-step probe for large grids, env-gated so normal test runs skip it.
    /// PROF_GRID="d,m" PROF_FOLDS=n cargo test --release -p signature-rs --lib bench_probe -- --nocapture
    /// Builds the series ONCE (the O(D^2) table build dominates), then times
    /// steady-state letter-displacement folds. RAYON_NUM_THREADS selects the
    /// regime under test.
    #[test]
    fn probe_large_grid_folds() {
        let Some(grid) = std::env::var_os("PROF_GRID") else { return };
        let grid = grid.to_str().unwrap().to_owned();
        let (d, m) = grid.split_once(',').unwrap();
        let (d, m) = (d.parse::<usize>().unwrap(), m.parse::<usize>().unwrap());
        let folds: usize = std::env::var("PROF_FOLDS").ok().and_then(|s| s.parse().ok()).unwrap_or(20);

        use ordered_float::NotNan;
        let builder = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(d)
            .with_max_degree(m);
        let t = std::time::Instant::now();
        let mut log_sig = builder.build::<NotNan<f64>>();
        let build = t.elapsed();
        println!("PROBE d={d} m={m} D={} series_build={build:?}", log_sig.series.coefficients.len());

        let basis = lyndon_rs::lyndon::LyndonBasis::<u8>::new(
            d,
            lyndon_rs::lyndon::Sort::Lexicographical,
        )
        .generate_basis(m);
        let mut seed = 0xfeed_u64;
        let mut lcg = |seed: &mut u64| {
            *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((*seed >> 33) as i64) as f64
        };
        let acc: Vec<NotNan<f64>> = (0..basis.len())
            .map(|_| NotNan::new(lcg(&mut seed)).unwrap())
            .collect();
        log_sig.series.coefficients.clone_from(&acc);
        let disp: Vec<NotNan<f64>> = (0..basis.len())
            .map(|k| NotNan::new(if k < d { lcg(&mut seed) } else { 0.0 }).unwrap())
            .collect();
        let mut seg = log_sig.clone();
        seg.series.coefficients.clone_from(&disp);

        // warm: the first fold takes the collecting-rebuild path, the rest
        // run the steady level-batch
        for _ in 0..4 {
            log_sig.concatenate_assign(&seg);
        }
        let t = std::time::Instant::now();
        for k in 0..folds {
            log_sig.concatenate_assign(&seg);
            std::hint::black_box(&log_sig);
            if k + 1 == folds {
                println!("PROBE d={d} m={m} threads={} fold={:?} ({} folds)",
                    std::thread::available_parallelism().map(|n| n.get()).unwrap_or(0),
                    t.elapsed() / (folds as u32), folds);
            }
        }
    }

}
