use super::*;

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
    ///
    /// # Caching
    ///
    /// The Lyndon basis, the structure-constant table, the BCH series and
    /// the commutator-DAG template are pure functions of the builder
    /// configuration (`T`, `U`, alphabet size, max degree), so they are
    /// built once per process and shared across every subsequent build:
    /// the first `build` for a configuration pays the (expensive) table
    /// construction, later builds are cheap clones of the cached plan.
    /// See [`LogSignature::concatenate_batch_coefficients`] for the fold
    /// machinery that consumes the plan.
    ///
    /// # Why `'static`
    ///
    /// The plan caches are process-wide and keyed by [`TypeId`], which
    /// requires `T: 'static` and `U: 'static`. All `Copy`/`Clone` numeric
    /// coefficient types (f32, f64, `Ratio`, `NotNan`, ...) qualify.
    ///
    /// [`TypeId`]: core::any::TypeId
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
    ) -> LogSignature<T, U>
    where
        T: 'static,
        U: 'static,
    {
        // The basis + structure-constant table are a pure function of the
        // builder configuration; build them once and share them across all
        // builds (see `series_plan_cache`). The BCH series + DAG template
        // are likewise a pure function of (max_degree, U) — the series is
        // cloned per build (Arc-shared internals, fresh coefficient vector)
        // and the DAG via `clone_shallow` (compiled structure + node-list
        // fixed point shared; fold scratch stays per instance).
        let plan = series_plan_cache::series_plan::<T, U>(&self.lyndon_basis, self.max_degree);
        let bch = series_plan_cache::bch_plan::<U>(self.max_degree);
        let basis_len = plan.basis.len();
        let coefficients = vec![U::default(); basis_len];
        let series = LieSeries::<T, U>::with_feasible_decompositions(
            Arc::clone(&plan.basis),
            coefficients,
            Arc::clone(&plan.table),
        );
        let bch_series = bch.series.clone();
        let dag = bch.dag.clone_shallow();
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
    ///
    /// For paths long enough to qualify for the batch driver's tournament
    /// reduction (see [`LogSignature::concatenate_batch_coefficients`]), the
    /// folds are reassociated along a balanced tree — a shape that depends
    /// only on the number of displacements, never on the thread count — so
    /// f32/f64 rounding may differ by a few ulps from strictly sequential
    /// folding, while exact coefficient types remain bit-identical.
    ///
    /// Like [`Self::build`], this caches the compiled plan (basis +
    /// structure-constant table + BCH/DAG template) process-wide per
    /// configuration, so repeated calls on one builder are cheap after the
    /// first; this is why `T` and `U` must be `'static` (the caches key on
    /// [`TypeId`]).
    ///
    /// [`TypeId`]: core::any::TypeId
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
            + Sub<Output = U>
            + 'static,
    >(
        &self,
        path: &ArrayView<U, D>,
    ) -> LogSignature<T, U>
    where
        T: 'static,
    {
        let mut log_sig = self.build();

        // One flat displacement buffer (each row is the elementwise
        // difference of two consecutive path points); the batch's
        // displacement slices borrow it. The per-window `Vec`s this
        // replaces (an owned difference array plus a `Vec` per segment)
        // were a measurable slice of the small-grid e2e's allocator
        // traffic.
        let rows = path.len_of(Axis(0));
        let n_displacements = rows.saturating_sub(1);
        let row_len = path.len().checked_div(rows).unwrap_or(0);
        let mut displacements = vec![U::default(); n_displacements * row_len];
        for (wi, window) in path.axis_windows(Axis(0), 2).into_iter().enumerate() {
            let start = window.index_axis(Axis(0), 0);
            let end = window.index_axis(Axis(0), 1);
            let dst = &mut displacements[wi * row_len..(wi + 1) * row_len];
            for (c, (e, s)) in dst.iter_mut().zip(end.iter().zip(start.iter())) {
                *c = e.clone() - s.clone();
            }
        }
        let slices: Vec<&[U]> = displacements
            .chunks(row_len.max(1))
            .take(n_displacements)
            .collect();
        log_sig.concatenate_batch_coefficients(&slices);

        log_sig
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::{Array2, array};
    use num_traits::ToPrimitive;
    use ordered_float::NotNan;
    use rstest::rstest;

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
}
