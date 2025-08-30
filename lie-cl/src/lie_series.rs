use std::marker::PhantomData;

use cubecl::{prelude::*, server::Handle, std::tensor::compact_strides};
use lie_rs::LieSeries;

pub struct LieSeriesCLData<R: Runtime, T: Numeric + CubeElement> {
    basis_map_shape: Vec<usize>,
    basis_coefficients_data: Handle,
    basis_map_data: Handle,
    basis_strides: Vec<usize>,
    coefficients_data: Handle,
    degrees_data: Handle,
    _t: PhantomData<T>,
    _r: PhantomData<R>,
}

impl<R: Runtime, T: Numeric + CubeElement> Clone for LieSeriesCLData<R, T> {
    fn clone(&self) -> Self {
        Self {
            basis_map_shape: self.basis_map_shape.clone(),
            basis_coefficients_data: self.basis_coefficients_data.clone(),
            basis_map_data: self.basis_map_data.clone(),
            basis_strides: self.basis_strides.clone(),
            coefficients_data: self.coefficients_data.clone(),
            degrees_data: self.degrees_data.clone(),
            _t: PhantomData,
            _r: PhantomData,
        }
    }
}

impl<R: Runtime, T: Numeric + CubeElement> LieSeriesCLData<R, T> {
    pub fn new(
        client: &ComputeClient<R::Server, R::Channel>,
        lie_series: &LieSeries<u8, impl Into<T> + Copy>,
    ) -> Self {
        let max_basis_terms = lie_series
            .commutator_basis_map
            .iter()
            .fold(0, |x, y| usize::max(x, y.len()));
        let basis_map_shape = vec![
            lie_series.commutator_basis.len(),
            lie_series.commutator_basis.len(),
            max_basis_terms,
        ];
        let mut basis_map_padded = vec![];
        let mut basis_coefficients_padded = vec![];

        for i in 0..lie_series.coefficients.len() {
            for j in 0..lie_series.coefficients.len() {
                let basis_indices =
                    &lie_series.commutator_basis_map[i * lie_series.basis.len() + j];
                let basis_coefficients =
                    &lie_series.commutator_basis_map_coefficients[i * lie_series.basis.len() + j];
                for k in 0..max_basis_terms {
                    if k >= basis_indices.len() {
                        basis_map_padded.push(-1i32);
                        basis_coefficients_padded.push(T::from_int(0));
                        continue;
                    }
                    basis_map_padded.push(basis_indices[k] as i32);
                    basis_coefficients_padded.push(basis_coefficients[k].into());
                }
            }
        }

        let basis_map_data = client.create(i32::as_bytes(&basis_map_padded));

        let basis_coefficients_data = client.create(T::as_bytes(&basis_coefficients_padded));

        let coefficients_data = client.create(T::as_bytes(
            &lie_series
                .coefficients
                .iter()
                .copied()
                .map(std::convert::Into::into)
                .collect::<Vec<_>>(),
        ));
        let degree_vec = lie_series
            .commutator_basis
            .iter()
            .map(|x| x.degree() as u32)
            .collect::<Vec<_>>();

        let degrees_data = client.create(u32::as_bytes(&degree_vec));
        let basis_strides = compact_strides(&basis_map_shape);

        Self {
            basis_map_shape,
            basis_map_data,
            basis_strides,
            basis_coefficients_data,
            coefficients_data,
            degrees_data,
            _t: PhantomData,
            _r: PhantomData,
        }
    }

    #[must_use]
    pub fn empty(
        client: &ComputeClient<R::Server, R::Channel>,
        lie_series: &LieSeries<u8, impl Into<T> + Copy>,
    ) -> Self {
        let max_basis_terms = lie_series
            .commutator_basis_map
            .iter()
            .fold(0, |x, y| usize::max(x, y.len()));
        let basis_map_shape = vec![
            lie_series.commutator_basis.len(),
            lie_series.commutator_basis.len(),
            max_basis_terms,
        ];
        let mut basis_map_padded = vec![];
        let mut basis_coefficients_padded = vec![];
        for i in 0..lie_series.coefficients.len() {
            for j in 0..lie_series.coefficients.len() {
                let basis_indices =
                    &lie_series.commutator_basis_map[i * lie_series.basis.len() + j];
                let basis_coefficients =
                    &lie_series.commutator_basis_map_coefficients[i * lie_series.basis.len() + j];
                for k in 0..max_basis_terms {
                    if k >= basis_indices.len() {
                        basis_map_padded.push(-1i32);
                        basis_coefficients_padded.push(T::from_int(0));
                        continue;
                    }
                    basis_map_padded.push(basis_indices[k] as i32);
                    basis_coefficients_padded.push(basis_coefficients[k].into());
                }
            }
        }
        let basis_map_data = client.create(i32::as_bytes(&basis_map_padded));

        let basis_coefficients_data = client.create(T::as_bytes(&basis_coefficients_padded));

        let coefficients_data = client.create(T::as_bytes(&vec![
            T::from_int(0);
            lie_series.coefficients.len()
        ]));
        let degree_vec = lie_series
            .commutator_basis
            .iter()
            .map(|x| x.degree() as u32)
            .collect::<Vec<_>>();

        let degrees_data = client.create(u32::as_bytes(&degree_vec));
        let basis_strides = compact_strides(&basis_map_shape);

        Self {
            basis_map_shape,
            basis_map_data,
            basis_strides,
            basis_coefficients_data,
            coefficients_data,
            degrees_data,
            _t: PhantomData,
            _r: PhantomData,
        }
    }

    #[must_use]
    pub fn kernel_arg<F: CubePrimitive>(&self) -> LieSeriesCLLaunch<F, R> {
        let basis_coefficients = unsafe {
            TensorArg::from_raw_parts::<F>(
                &self.basis_coefficients_data,
                &self.basis_strides,
                &self.basis_map_shape,
                1,
            )
        };

        let basis_map = unsafe {
            TensorArg::from_raw_parts::<i32>(
                &self.basis_map_data,
                &self.basis_strides,
                &self.basis_map_shape,
                1,
            )
        };

        let coefficients = unsafe {
            ArrayArg::from_raw_parts::<T>(&self.coefficients_data, self.basis_map_shape[0], 1)
        };

        let degrees = unsafe {
            ArrayArg::from_raw_parts::<u32>(&self.degrees_data, self.basis_map_shape[0], 1)
        };
        LieSeriesCLLaunch::new(basis_coefficients, basis_map, coefficients, degrees)
    }

    pub fn read_coefficients(self, client: &ComputeClient<R::Server, R::Channel>) -> Vec<T> {
        let bytes = client.read_one(self.coefficients_data.binding());
        T::from_bytes(&bytes).to_vec()
    }
}

#[derive(CubeType, CubeLaunch)]
pub struct LieSeriesCL<T: CubePrimitive> {
    /// An N * N * C matrix where the first two dimensions represent the
    /// i, j pair and the last dimension represents the basis coefficient
    /// for each index in the basis map
    pub basis_coefficients: Tensor<T>,
    /// The N * N * C matrix of where to locate the basis coefficients
    /// for a given i, j pair. Values of `-1` mean padding and should
    /// be ignored.
    pub basis_map: Tensor<i32>,
    /// The coefficients of the Lie Series
    pub coefficients: Array<T>,
    /// The degrees for each term in the Lie Series
    pub degrees: Array<u32>,
}

#[cube(launch_unchecked)]
pub fn lie_series_commutation<F: Float>(
    input_a: &LieSeriesCL<F>,
    input_b: &LieSeriesCL<F>,
    output: &mut LieSeriesCL<Atomic<F>>,
) {
    let max_degree = input_a.degrees[input_a.degrees.len() - 1];
    for i in 0..input_a.basis_map.shape(0) {
        for j in 0..input_a.basis_map.shape(1) {
            let coeff = if i == j || input_a.degrees[i] + input_b.degrees[j] > max_degree {
                F::new(0.0)
            } else {
                input_a.coefficients[i] * input_b.coefficients[j]
            };

            for k in 0..input_a.basis_map.shape(2) {
                let basis_index = input_a.basis_map
                    [i * input_a.basis_map.stride(0) + j * input_a.basis_map.stride(1) + k];
                let basis_coefficient = input_a.basis_coefficients
                    [i * input_a.basis_map.stride(0) + j * input_a.basis_map.stride(1) + k];
                if basis_index > -1 {
                    Atomic::add(
                        &output.coefficients[basis_index as u32],
                        basis_coefficient * coeff,
                    );
                }
            }
        }
    }
}

#[cube(launch_unchecked)]
pub fn lie_series_commutation_3d<F: Float>(
    input_a: &LieSeriesCL<F>,
    input_b: &LieSeriesCL<F>,
    output: &mut LieSeriesCL<Atomic<F>>,
) {
    let i = CUBE_POS_X * CUBE_DIM_X + UNIT_POS_X;
    let j = CUBE_POS_Y * CUBE_DIM_Y + UNIT_POS_Y;
    let k = CUBE_POS_Z * CUBE_DIM_Z + UNIT_POS_Z;

    if i >= input_a.basis_map.shape(0)
        || j >= input_a.basis_map.shape(1)
        || k >= input_a.basis_map.shape(2)
    {
        terminate!();
    }

    let max_degree = input_a.degrees[input_a.degrees.len() - 1];
    if i == j || input_a.degrees[i] + input_b.degrees[j] > max_degree {
        terminate!();
    }

    let basis_index =
        input_a.basis_map[i * input_a.basis_map.stride(0) + j * input_a.basis_map.stride(1) + k];
    if basis_index == -1 {
        terminate!();
    }

    let coeff = input_a.coefficients[i] * input_b.coefficients[j];
    if coeff == F::from_int(0) {
        terminate!();
    }
    let basis_coefficient = input_a.basis_coefficients
        [i * input_a.basis_map.stride(0) + j * input_a.basis_map.stride(1) + k];
    if basis_coefficient == F::from_int(0) {
        terminate!();
    }

    Atomic::add(
        &output.coefficients[basis_index as u32],
        basis_coefficient * coeff,
    );
}

#[cfg(test)]
mod test {
    use std::cmp::min;

    use cubecl::{
        benchmark::Benchmark,
        cuda::CudaRuntime,
        future,
        wgpu::{WgpuDevice, WgpuRuntime},
    };
    use lyndon_rs::{LyndonBasis, Sort};
    use ordered_float::NotNan;
    use rstest::rstest;

    use super::*;

    pub struct SerialCommutationBench<R: Runtime> {
        pub num_generators: usize,
        pub basis_depth: usize,
        pub a_coefficients: Vec<f32>,
        pub b_coefficients: Vec<f32>,
        pub lie_series: LieSeries<u8, NotNan<f32>>,
        client: ComputeClient<R::Server, R::Channel>,
    }

    impl<R: Runtime> Benchmark for SerialCommutationBench<R> {
        type Input = (LieSeriesCLData<R, f32>, LieSeriesCLData<R, f32>);
        type Output = LieSeriesCLData<R, f32>;

        fn prepare(&self) -> Self::Input {
            let basis = LyndonBasis::<u8>::new(self.num_generators, Sort::Lexicographical)
                .generate_basis(self.basis_depth);

            let a = LieSeries::new(
                basis.clone(),
                self.a_coefficients
                    .iter()
                    .copied()
                    .map(|x| NotNan::try_from(x).unwrap())
                    .collect::<Vec<_>>(),
            );
            let b = LieSeries::new(
                basis.clone(),
                self.b_coefficients
                    .iter()
                    .copied()
                    .map(|x| NotNan::try_from(x).unwrap())
                    .collect::<Vec<_>>(),
            );
            let input_a = LieSeriesCLData::<R, f32>::new(&self.client, &a);
            let input_b = LieSeriesCLData::<R, f32>::new(&self.client, &b);
            (input_a, input_b)
        }

        fn name(&self) -> String {
            format!(
                "{}-commutation-{:?}-{:?}",
                R::name(&self.client),
                self.num_generators,
                self.basis_depth
            )
            .to_lowercase()
        }

        fn sync(&self) {
            future::block_on(self.client.sync());
        }

        fn execute(&self, input: Self::Input) -> Result<Self::Output, String> {
            let output = LieSeriesCLData::<R, f32>::empty(&self.client, &self.lie_series);

            unsafe {
                lie_series_commutation::launch_unchecked(
                    &self.client,
                    CubeCount::Static(1, 1, 1),
                    CubeDim::new(1, 1, 1),
                    input.0.kernel_arg::<f32>(),
                    input.1.kernel_arg(),
                    output.kernel_arg(),
                );
            }
            Ok(output)
        }
    }

    pub struct ParallelCommutationBench<R: Runtime> {
        pub num_generators: usize,
        pub basis_depth: usize,
        pub a_coefficients: Vec<f32>,
        pub b_coefficients: Vec<f32>,
        pub lie_series: LieSeries<u8, NotNan<f32>>,
        client: ComputeClient<R::Server, R::Channel>,
    }

    impl<R: Runtime> Benchmark for ParallelCommutationBench<R> {
        type Input = (LieSeriesCLData<R, f32>, LieSeriesCLData<R, f32>);
        type Output = LieSeriesCLData<R, f32>;

        fn prepare(&self) -> Self::Input {
            let basis = LyndonBasis::<u8>::new(self.num_generators, Sort::Lexicographical)
                .generate_basis(self.basis_depth);

            let a = LieSeries::new(
                basis.clone(),
                self.a_coefficients
                    .iter()
                    .copied()
                    .map(|x| NotNan::try_from(x).unwrap())
                    .collect::<Vec<_>>(),
            );
            let b = LieSeries::new(
                basis.clone(),
                self.b_coefficients
                    .iter()
                    .copied()
                    .map(|x| NotNan::try_from(x).unwrap())
                    .collect::<Vec<_>>(),
            );
            let input_a = LieSeriesCLData::<R, f32>::new(&self.client, &a);
            let input_b = LieSeriesCLData::<R, f32>::new(&self.client, &b);
            (input_a, input_b)
        }

        fn name(&self) -> String {
            format!(
                "{}-commutation-{:?}-{:?}",
                R::name(&self.client),
                self.num_generators,
                self.basis_depth
            )
            .to_lowercase()
        }

        fn sync(&self) {
            future::block_on(self.client.sync());
        }

        fn execute(&self, input: Self::Input) -> Result<Self::Output, String> {
            let output = LieSeriesCLData::<R, f32>::empty(&self.client, &self.lie_series);
            let dim_i = input.0.basis_map_shape[0] as u32;
            let dim_j = input.0.basis_map_shape[1] as u32;
            let dim_k = input.0.basis_map_shape[2] as u32;

            let cube_x = min(dim_i, 32);
            let cube_y = min(dim_j, 32);
            let cube_z = min(dim_k, 1024 / (cube_x * cube_y));

            let cubes_x = dim_i.div_ceil(cube_x);
            let cubes_y = dim_j.div_ceil(cube_y);
            let cubes_z = dim_k.div_ceil(cube_z);

            let cube_count = CubeCount::Static(cubes_x, cubes_y, cubes_z);
            let cube_dim = CubeDim::new(cube_x, cube_y, cube_z);

            unsafe {
                lie_series_commutation_3d::launch_unchecked(
                    &self.client,
                    cube_count,
                    cube_dim,
                    input.0.kernel_arg::<f32>(),
                    input.1.kernel_arg(),
                    output.kernel_arg(),
                );
            }
            Ok(output)
        }
    }

    #[rstest]
    #[case(2, 2, vec![1., 2., 3.], vec![4., 5., 6.], vec![0., 0., -3.])]
    #[case(2, 2, vec![3., 2., 1.], vec![1., 2., 3.], vec![0., 0., 4.])]
    #[case(2, 3, vec![1., 2., 3., 4., 5.], vec![6., 7., 8., 9., 10.], vec![0., 0., -5., -10., 5.])]
    #[case(3, 3,
        vec![1., 2., 3., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
        vec![5., 3., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
        vec![0., 0., 0., -7., -14., -7., 0., 0., 0., 0., 0., 0., 0., 0.])]
    #[case(3, 3,
        vec![1., 2., 3., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
        vec![0., 0., 0., -7., -14., -7., 0., 0., 0., 0., 0., 0., 0., 0.],
        vec![0., 0., 0., 0., 0., 0., -7., -14., 14., 14., 49., 42., -14., 21.])]
    #[case(3, 4,
        vec![
            1., 2., 3., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
        ],
        vec![
            5., 3., 1., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
        ],
        vec![
            0., 0., 0., -7., -14., -7., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
        ],
    )]
    #[case(3, 4, vec![
            1., 2., 3., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
        ],
        vec![
            0., 0., 0., -7., -14., -7., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
        ],
        vec![
            0., 0., 0., 0., 0., 0., -7., -14.,
            14., 14., 49., 42., -14., 21., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
        ])]
    fn test_lie_series_commutation_on_gpu(
        #[case] num_generators: usize,
        #[case] basis_depth: usize,
        #[case] a_coefficients: Vec<f32>,
        #[case] b_coefficients: Vec<f32>,
        #[case] expected_coefficients: Vec<f32>,
    ) {
        let client = WgpuRuntime::client(&WgpuDevice::default());
        let basis = LyndonBasis::<u8>::new(num_generators, Sort::Lexicographical)
            .generate_basis(basis_depth);
        let a_coefficients = a_coefficients
            .into_iter()
            .map(|x| NotNan::try_from(x).unwrap())
            .collect::<Vec<_>>();
        let b_coefficients = b_coefficients
            .into_iter()
            .map(|x| NotNan::try_from(x).unwrap())
            .collect::<Vec<_>>();
        let a = LieSeries::new(basis.clone(), a_coefficients);
        let b = LieSeries::new(basis, b_coefficients);
        let input_a = LieSeriesCLData::<WgpuRuntime, f32>::new(&client, &a);
        let input_b = LieSeriesCLData::<WgpuRuntime, f32>::new(&client, &b);
        let output = LieSeriesCLData::<WgpuRuntime, f32>::empty(&client, &a);

        unsafe {
            lie_series_commutation::launch_unchecked(
                &client,
                CubeCount::Static(1, 1, 1),
                CubeDim::new(1, 1, 1),
                input_a.kernel_arg::<f32>(),
                input_b.kernel_arg(),
                output.kernel_arg(),
            );
        }

        let coefficients = output.read_coefficients(&client);
        assert_eq!(coefficients.len(), expected_coefficients.len());
        dbg!(&coefficients);
        for (c, e_c) in coefficients
            .into_iter()
            .zip(expected_coefficients.into_iter())
        {
            assert!((c - e_c).abs() < 0.001, "{c} != {e_c}");
        }
    }

    #[rstest]
    #[case(2, 2, vec![1., 2., 3.], vec![4., 5., 6.], vec![0., 0., -3.])]
    #[case(2, 2, vec![3., 2., 1.], vec![1., 2., 3.], vec![0., 0., 4.])]
    #[case(2, 3, vec![1., 2., 3., 4., 5.], vec![6., 7., 8., 9., 10.], vec![0., 0., -5., -10., 5.])]
    #[case(3, 3,
        vec![1., 2., 3., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
        vec![5., 3., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
        vec![0., 0., 0., -7., -14., -7., 0., 0., 0., 0., 0., 0., 0., 0.])]
    #[case(3, 3,
        vec![1., 2., 3., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
        vec![0., 0., 0., -7., -14., -7., 0., 0., 0., 0., 0., 0., 0., 0.],
        vec![0., 0., 0., 0., 0., 0., -7., -14., 14., 14., 49., 42., -14., 21.])]
    #[case(3, 4,
        vec![
            1., 2., 3., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
        ],
        vec![
            5., 3., 1., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
        ],
        vec![
            0., 0., 0., -7., -14., -7., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
        ],
    )]
    #[case(3, 4, vec![
            1., 2., 3., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
        ],
        vec![
            0., 0., 0., -7., -14., -7., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
        ],
        vec![
            0., 0., 0., 0., 0., 0., -7., -14.,
            14., 14., 49., 42., -14., 21., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0.,
        ])]
    fn test_lie_series_commutation_3d(
        #[case] num_generators: usize,
        #[case] basis_depth: usize,
        #[case] a_coefficients: Vec<f32>,
        #[case] b_coefficients: Vec<f32>,
        #[case] expected_coefficients: Vec<f32>,
    ) {
        use std::cmp::min;

        let client = WgpuRuntime::client(&WgpuDevice::default());
        let basis = LyndonBasis::<u8>::new(num_generators, Sort::Lexicographical)
            .generate_basis(basis_depth);
        let a_coefficients = a_coefficients
            .into_iter()
            .map(|x| NotNan::try_from(x).unwrap())
            .collect::<Vec<_>>();
        let b_coefficients = b_coefficients
            .into_iter()
            .map(|x| NotNan::try_from(x).unwrap())
            .collect::<Vec<_>>();
        let a = LieSeries::new(basis.clone(), a_coefficients);
        let b = LieSeries::new(basis, b_coefficients);
        let input_a = LieSeriesCLData::<WgpuRuntime, f32>::new(&client, &a);
        let input_b = LieSeriesCLData::<WgpuRuntime, f32>::new(&client, &b);
        let output = LieSeriesCLData::<WgpuRuntime, f32>::empty(&client, &a);

        let dim_i = input_a.basis_map_shape[0] as u32;
        let dim_j = input_a.basis_map_shape[1] as u32;
        let dim_k = input_a.basis_map_shape[2] as u32;

        let cube_x = min(dim_i, 32);
        let cube_y = min(dim_j, 32);
        let cube_z = min(dim_k, 1024 / (cube_x * cube_y));

        let cubes_x = dim_i.div_ceil(cube_x);
        let cubes_y = dim_j.div_ceil(cube_y);
        let cubes_z = dim_k.div_ceil(cube_z);

        let cube_count = CubeCount::Static(cubes_x, cubes_y, cubes_z);
        let cube_dim = CubeDim::new(cube_x, cube_y, cube_z);

        unsafe {
            lie_series_commutation_3d::launch_unchecked(
                &client,
                cube_count,
                cube_dim,
                input_a.kernel_arg::<f32>(),
                input_b.kernel_arg(),
                output.kernel_arg(),
            );
        }

        let coefficients = output.read_coefficients(&client);
        assert_eq!(coefficients.len(), expected_coefficients.len());
        dbg!(&coefficients);
        for (c, e_c) in coefficients
            .into_iter()
            .zip(expected_coefficients.into_iter())
        {
            assert!((c - e_c).abs() < 0.001, "{c} != {e_c}");
        }
    }

    #[rstest]
    #[case(3, 5)]
    #[case(4, 5)]
    #[case(5, 5)]
    #[case(6, 4)]
    #[case(7, 4)]
    #[case(8, 4)]
    fn test_serial_gpu_benchmark(#[case] num_generators: usize, #[case] basis_depth: usize) {
        let basis = LyndonBasis::<u8>::new(num_generators, Sort::Lexicographical)
            .generate_basis(basis_depth);
        let a_coefficients = (0..basis.len())
            .map(|x| NotNan::<f32>::try_from(x as f32).unwrap())
            .collect::<Vec<_>>();
        let b_coefficients = (0..basis.len()).rev().map(|x| x as f32).collect::<Vec<_>>();
        let client = WgpuRuntime::client(&WgpuDevice::default());
        let lie_series = LieSeries::new(basis, a_coefficients.clone());
        let a_coefficients = a_coefficients.into_iter().map(|x| x.into_inner()).collect();
        let bench = SerialCommutationBench::<WgpuRuntime> {
            num_generators,
            basis_depth,
            a_coefficients,
            b_coefficients,
            lie_series,
            client,
        };

        println!("{}", bench.name());
        println!(
            "{}",
            bench.run(cubecl::benchmark::TimingMethod::System).unwrap()
        );
    }

    #[rstest]
    #[case(3, 5)]
    #[case(4, 5)]
    #[case(5, 5)]
    #[case(6, 4)]
    #[case(7, 4)]
    #[case(8, 4)]
    fn test_parallel_gpu_benchmark(#[case] num_generators: usize, #[case] basis_depth: usize) {
        use cubecl::cuda::CudaDevice;

        let basis = LyndonBasis::<u8>::new(num_generators, Sort::Lexicographical)
            .generate_basis(basis_depth);
        let a_coefficients = (0..basis.len())
            .map(|x| NotNan::<f32>::try_from(x as f32).unwrap())
            .collect::<Vec<_>>();
        let b_coefficients = (0..basis.len()).rev().map(|x| x as f32).collect::<Vec<_>>();
        let client = CudaRuntime::client(&CudaDevice::default());
        let lie_series = LieSeries::new(basis, a_coefficients.clone());
        let a_coefficients = a_coefficients.into_iter().map(|x| x.into_inner()).collect();
        let bench = ParallelCommutationBench::<CudaRuntime> {
            num_generators,
            basis_depth,
            a_coefficients,
            b_coefficients,
            lie_series,
            client,
        };

        println!("{}", bench.name());
        println!(
            "{}",
            bench.run(cubecl::benchmark::TimingMethod::System).unwrap()
        );
    }
}
