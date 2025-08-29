use std::fmt::Debug;
use std::marker::PhantomData;

use cubecl::{prelude::*, server::Handle, std::tensor::compact_strides};
use lie_rs::LieSeries;

pub struct LieSeriesCLData<R: Runtime, F: Float + CubeElement> {
    basis_map_shape: Vec<usize>,
    basis_coefficients_data: Handle,
    basis_map_data: Handle,
    basis_strides: Vec<usize>,
    coefficients_data: Handle,
    degrees_data: Handle,
    _f: PhantomData<F>,
    _r: PhantomData<R>,
}

impl<R: Runtime, F: Float + CubeElement> Clone for LieSeriesCLData<R, F> {
    fn clone(&self) -> Self {
        Self {
            basis_map_shape: self.basis_map_shape.clone(),
            basis_coefficients_data: self.basis_coefficients_data.clone(),
            basis_map_data: self.basis_map_data.clone(),
            basis_strides: self.basis_strides.clone(),
            coefficients_data: self.coefficients_data.clone(),
            degrees_data: self.degrees_data.clone(),
            _f: PhantomData,
            _r: PhantomData,
        }
    }
}

impl<R: Runtime, F: Float + CubeElement> LieSeriesCLData<R, F> {
    pub fn new(
        client: &ComputeClient<R::Server, R::Channel>,
        lie_series: &LieSeries<u8, impl Into<F> + Copy>,
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
                        basis_coefficients_padded.push(F::new(0.0));
                        continue;
                    }
                    basis_map_padded.push(basis_indices[k] as i32);
                    basis_coefficients_padded.push(basis_coefficients[k].into());
                }
            }
        }

        let basis_map_data = client.create(i32::as_bytes(&basis_map_padded));

        let basis_coefficients_data = client.create(F::as_bytes(&basis_coefficients_padded));

        let coefficients_data = client.create(F::as_bytes(
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
            _f: PhantomData,
            _r: PhantomData,
        }
    }

    #[must_use]
    pub fn empty(
        client: &ComputeClient<R::Server, R::Channel>,
        lie_series: &LieSeries<u8, impl Into<F> + Copy>,
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
                        basis_coefficients_padded.push(F::new(0.0));
                        continue;
                    }
                    basis_map_padded.push(basis_indices[k] as i32);
                    basis_coefficients_padded.push(basis_coefficients[k].into());
                }
            }
        }
        let basis_map_data = client.create(i32::as_bytes(&basis_map_padded));

        let basis_coefficients_data = client.create(F::as_bytes(&basis_coefficients_padded));

        let coefficients_data = client.create(F::as_bytes(&vec![
            F::new(0.0);
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
            _f: PhantomData,
            _r: PhantomData,
        }
    }

    #[must_use]
    pub fn kernel_arg(&self) -> LieSeriesCLLaunch<F, R> {
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
            ArrayArg::from_raw_parts::<F>(&self.coefficients_data, self.basis_map_shape[0], 1)
        };

        let degrees = unsafe {
            ArrayArg::from_raw_parts::<u32>(&self.degrees_data, self.basis_map_shape[0], 1)
        };
        LieSeriesCLLaunch::new(basis_coefficients, basis_map, coefficients, degrees)
    }

    pub fn read_coefficients(self, client: &ComputeClient<R::Server, R::Channel>) -> Vec<F> {
        let bytes = client.read_one(self.coefficients_data.binding());
        F::from_bytes(&bytes).to_vec()
    }
}

#[derive(CubeType, CubeLaunch)]
pub struct LieSeriesCL<F: Float> {
    /// An N * N * C matrix where the first two dimensions represent the
    /// i, j pair and the last dimension represents the basis coefficient
    /// for each index in the basis map
    pub basis_coefficients: Tensor<F>,
    /// The N * N * C matrix of where to locate the basis coefficients
    /// for a given i, j pair. Values of `-1` mean padding and should
    /// be ignored.
    pub basis_map: Tensor<i32>,
    /// The coefficients of the Lie Series
    pub coefficients: Array<F>,
    /// The degrees for each term in the Lie Series
    pub degrees: Array<u32>,
}

#[cube(launch_unchecked)]
pub fn lie_series_commutation<F: Float>(
    input_a: &LieSeriesCL<F>,
    input_b: &LieSeriesCL<F>,
    output: &mut LieSeriesCL<F>,
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
                    output.coefficients[basis_index as u32] += basis_coefficient * coeff;
                }
            }
        }
    }
}

#[cfg(test)]
mod test {
    use cubecl::wgpu::{WgpuDevice, WgpuRuntime};
    use lyndon_rs::{LyndonBasis, Sort};
    use ordered_float::NotNan;
    use rstest::rstest;

    use super::*;

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
}
