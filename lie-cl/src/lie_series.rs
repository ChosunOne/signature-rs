use std::marker::PhantomData;

use cubecl::{prelude::*, server::Handle, std::tensor::compact_strides, zspace::Strides};
use lie_rs::LieSeries;

const VECTORIZATION: usize = 4;

pub struct LieSeriesCLData<R, T> {
    pub basis_map_shape: Vec<usize>,
    basis_coefficients_data: Handle,
    basis_map_data: Handle,
    basis_strides: Strides,
    pub coefficients_data: Handle,
    coefficients_len: usize,
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
            coefficients_len: self.coefficients_len,
            degrees_data: self.degrees_data.clone(),
            _t: PhantomData,
            _r: PhantomData,
        }
    }
}

impl<R: Runtime, T: Numeric + CubeElement> LieSeriesCLData<R, T> {
    pub fn new(client: &ComputeClient<R>, lie_series: &LieSeries<u8, impl Into<T> + Copy>) -> Self {
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

        let basis_map_data = client.create_from_slice(i32::as_bytes(&basis_map_padded));

        let basis_coefficients_data =
            client.create_from_slice(T::as_bytes(&basis_coefficients_padded));

        let coefficients = lie_series
            .coefficients
            .iter()
            .copied()
            .map(Into::into)
            .collect::<Vec<T>>();

        let coefficients_len = coefficients.len().div_ceil(VECTORIZATION) * VECTORIZATION;

        let mut coefficients = coefficients;
        coefficients.resize(coefficients_len, T::from_int(0));
        let coefficients_data = client.create_from_slice(T::as_bytes(&coefficients));

        let degree_vec = lie_series
            .commutator_basis
            .iter()
            .map(|x| x.degree() as u32)
            .collect::<Vec<_>>();

        let degrees_data = client.create_from_slice(u32::as_bytes(&degree_vec));
        let basis_strides = compact_strides(&basis_map_shape);

        Self {
            basis_map_shape,
            basis_map_data,
            basis_strides,
            basis_coefficients_data,
            coefficients_data,
            coefficients_len,
            degrees_data,
            _t: PhantomData,
            _r: PhantomData,
        }
    }

    #[must_use]
    pub fn empty(
        client: &ComputeClient<R>,
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
        let basis_map_data = client.create_from_slice(i32::as_bytes(&basis_map_padded));

        let basis_coefficients_data =
            client.create_from_slice(T::as_bytes(&basis_coefficients_padded));

        let coefficients_len =
            lie_series.coefficients.len().div_ceil(VECTORIZATION) * VECTORIZATION;

        let coefficients_data =
            client.create_from_slice(T::as_bytes(&vec![T::from_int(0); coefficients_len]));

        let degree_vec = lie_series
            .commutator_basis
            .iter()
            .map(|x| x.degree() as u32)
            .collect::<Vec<_>>();

        let degrees_data = client.create_from_slice(u32::as_bytes(&degree_vec));
        let basis_strides = compact_strides(&basis_map_shape);

        Self {
            basis_map_shape,
            basis_map_data,
            basis_strides,
            basis_coefficients_data,
            coefficients_data,
            coefficients_len,
            degrees_data,
            _t: PhantomData,
            _r: PhantomData,
        }
    }

    #[must_use]
    pub fn kernel_arg<E: CubePrimitive>(&self) -> LieSeriesCLLaunch<E, R> {
        let basis_coefficients = unsafe {
            TensorArg::from_raw_parts(
                self.basis_coefficients_data.clone(),
                self.basis_strides.clone(),
                self.basis_map_shape.clone().into(),
            )
        };

        let basis_map = unsafe {
            TensorArg::from_raw_parts(
                self.basis_map_data.clone(),
                self.basis_strides.clone(),
                self.basis_map_shape.clone().into(),
            )
        };

        let coefficients = unsafe {
            ArrayArg::from_raw_parts(self.coefficients_data.clone(), self.coefficients_len)
        };

        let degrees =
            unsafe { ArrayArg::from_raw_parts(self.degrees_data.clone(), self.basis_map_shape[0]) };
        LieSeriesCLLaunch::new(basis_coefficients, basis_map, coefficients, degrees)
    }

    pub fn read_coefficients(self, client: &ComputeClient<R>) -> Vec<T> {
        let n = self.basis_map_shape[0];
        let bytes = client.read_one(self.coefficients_data).unwrap();
        T::from_bytes(&bytes).iter().take(n).copied().collect()
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
pub fn lie_series_scalar_multiplication<F: Numeric + CubeElement, N: Size>(
    input_a: &LieSeriesCL<Vector<F, N>>,
    scalar: F,
    output: &mut LieSeriesCL<Vector<F, N>>,
) {
    let i = ABSOLUTE_POS;
    if i < input_a.coefficients.len() {
        let scalar = Vector::new(scalar);
        output.coefficients[i] = input_a.coefficients[i] * scalar;
    }
}

#[cube(launch_unchecked)]
pub fn lie_series_addition<F: Numeric, N: Size>(
    input_a: &LieSeriesCL<Vector<F, N>>,
    input_b: &LieSeriesCL<Vector<F, N>>,
    output: &mut LieSeriesCL<Vector<F, N>>,
) {
    let i = ABSOLUTE_POS;
    if i < input_a.coefficients.len() {
        output.coefficients[i] = input_a.coefficients[i] + input_b.coefficients[i];
    }
}

#[cube(launch_unchecked)]
pub fn lie_series_commutation<F: Numeric>(
    input_a: &LieSeriesCL<F>,
    input_b: &LieSeriesCL<F>,
    output: &mut LieSeriesCL<Atomic<F>>,
) {
    let i = (CUBE_POS_X * CUBE_DIM_X + UNIT_POS_X) as usize;
    let j = (CUBE_POS_Y * CUBE_DIM_Y + UNIT_POS_Y) as usize;
    let k = (CUBE_POS_Z * CUBE_DIM_Z + UNIT_POS_Z) as usize;

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

    output.coefficients[basis_index as usize].fetch_add(basis_coefficient * coeff);
}

#[cfg(test)]
mod test {
    use std::cmp::min;

    use cubecl::{
        benchmark::{Benchmark, TimingMethod},
        cuda::CudaRuntime,
        future,
        wgpu::{WgpuDevice, WgpuRuntime},
    };
    use lyndon_rs::{LyndonBasis, Sort};
    use ordered_float::NotNan;
    use rstest::rstest;

    use super::*;

    pub struct ParallelCommutationBench<R: Runtime> {
        pub num_generators: usize,
        pub basis_depth: usize,
        pub a_coefficients: Vec<f32>,
        pub b_coefficients: Vec<f32>,
        pub lie_series: LieSeries<u8, NotNan<f32>>,
        client: ComputeClient<R>,
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
            let _ = future::block_on(self.client.sync());
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
            let cube_dim = CubeDim::new_3d(cube_x, cube_y, cube_z);

            unsafe {
                lie_series_commutation::launch_unchecked::<f32, R>(
                    &self.client,
                    cube_count,
                    cube_dim,
                    input.0.kernel_arg(),
                    input.1.kernel_arg(),
                    output.kernel_arg(),
                );
            }
            Ok(output)
        }
    }

    #[rstest]
    fn test_scalar_multiplication() {
        use cubecl::cuda::CudaDevice;

        let basis = LyndonBasis::<u8>::new(3, Sort::Lexicographical).generate_basis(6);
        let a_coefficients = (0..basis.len())
            .map(|x| NotNan::<f32>::try_from(x as f32).unwrap())
            .collect::<Vec<_>>();
        let lie_series = LieSeries::new(basis, a_coefficients.clone());
        let client = CudaRuntime::client(&CudaDevice::default());
        let lie_series_data = LieSeriesCLData::<CudaRuntime, f32>::new(&client, &lie_series);
        let line_size = 4;

        unsafe {
            lie_series_scalar_multiplication::launch_unchecked::<f32, CudaRuntime>(
                &client,
                CubeCount::Static(1, 1, 1),
                CubeDim::new_3d(a_coefficients.len() as u32, 1, 1),
                line_size as usize,
                lie_series_data.kernel_arg(),
                2.0f32,
                lie_series_data.kernel_arg(),
            );
        }

        let coefficients = lie_series_data.read_coefficients(&client);
        assert_eq!(coefficients.len(), a_coefficients.len());

        for (c, o_c) in coefficients.into_iter().zip(
            a_coefficients
                .into_iter()
                .map(ordered_float::NotNan::into_inner),
        ) {
            assert!((c - 2. * o_c).abs() < 0.001, "{c} != 2.0 * {o_c}");
        }
    }

    #[rstest]
    fn test_addition() {
        use cubecl::cuda::CudaDevice;

        let basis = LyndonBasis::<u8>::new(3, Sort::Lexicographical).generate_basis(6);
        let a_coefficients = (0..basis.len())
            .map(|x| NotNan::<f32>::try_from(x as f32).unwrap())
            .collect::<Vec<_>>();
        let lie_series = LieSeries::new(basis, a_coefficients.clone());
        let client = CudaRuntime::client(&CudaDevice::default());
        let lie_series_data = LieSeriesCLData::<CudaRuntime, f32>::new(&client, &lie_series);
        let line_size = 4;
        let output = LieSeriesCLData::<CudaRuntime, f32>::empty(&client, &lie_series);

        unsafe {
            lie_series_addition::launch_unchecked::<f32, CudaRuntime>(
                &client,
                CubeCount::Static(1, 1, 1),
                CubeDim::new_3d(a_coefficients.len() as u32, 1, 1),
                line_size as usize,
                lie_series_data.kernel_arg(),
                lie_series_data.kernel_arg(),
                output.kernel_arg(),
            );
        }

        let coefficients = output.read_coefficients(&client);
        assert_eq!(coefficients.len(), a_coefficients.len());

        for (c, o_c) in coefficients.into_iter().zip(
            a_coefficients
                .into_iter()
                .map(ordered_float::NotNan::into_inner),
        ) {
            assert!((c - 2. * o_c).abs() < 0.001, "{c} != 2.0 * {o_c}");
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
    fn test_lie_series_commutation(
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
        let cube_dim = CubeDim::new_3d(cube_x, cube_y, cube_z);

        unsafe {
            lie_series_commutation::launch_unchecked::<f32, WgpuRuntime>(
                &client,
                cube_count,
                cube_dim,
                input_a.kernel_arg(),
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
        let a_coefficients = a_coefficients
            .into_iter()
            .map(ordered_float::NotNan::into_inner)
            .collect();
        let bench = ParallelCommutationBench::<CudaRuntime> {
            num_generators,
            basis_depth,
            a_coefficients,
            b_coefficients,
            lie_series,
            client,
        };

        println!("{}", bench.name());
        println!("{}", bench.run(TimingMethod::System).unwrap());
    }
}
