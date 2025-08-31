use std::collections::HashMap;

use cubecl::{prelude::*, server::Handle, std::tensor::compact_strides};
use lie_cl::lie_series::{LieSeriesCL, LieSeriesCLData};
use signature_rs::{CommutatorTerm, LogSignature};

pub struct LogSignatureCLData<R: Runtime, T: Numeric + CubeElement> {
    /// The necessary information for calculating Lie Series commutations
    lie_series_data: LieSeriesCLData<R, T>,
    /// The coefficients for each of the terms in the BCH series
    bch_coefficients_data: Handle,
    /// The number of coefficients in the `bch_coefficients_data`
    bch_coefficients_len: usize,
    /// The dependencies for each term in the series, e.g.
    /// A = [0, -1], B = [1, -1], [A, B] = [0, 1], [A, [A, B]] = [0, 2], ...
    bch_commutator_basis_data: Handle,
    bch_commutator_basis_shape: Vec<usize>,
    bch_commutator_basis_stride: Vec<usize>,
}

impl<R: Runtime, T: Numeric + CubeElement> Clone for LogSignatureCLData<R, T> {
    fn clone(&self) -> Self {
        Self {
            lie_series_data: self.lie_series_data.clone(),
            bch_coefficients_data: self.bch_coefficients_data.clone(),
            bch_coefficients_len: self.bch_coefficients_len,
            bch_commutator_basis_data: self.bch_commutator_basis_data.clone(),
            bch_commutator_basis_shape: self.bch_commutator_basis_shape.clone(),
            bch_commutator_basis_stride: self.bch_commutator_basis_stride.clone(),
        }
    }
}

impl<R: Runtime, T: Numeric + CubeElement> LogSignatureCLData<R, T> {
    pub fn new(
        client: &ComputeClient<R::Server, R::Channel>,
        log_signature: &LogSignature<u8, T>,
    ) -> Self {
        let lie_series_data = LieSeriesCLData::<R, T>::new(client, &log_signature.series);
        let bch_coefficients_data =
            client.create(T::as_bytes(&log_signature.bch_series.coefficients));

        let mut bch_commutator_basis =
            Vec::with_capacity(log_signature.bch_series.commutator_basis.len());
        let bch_commutator_basis_shape = vec![log_signature.bch_series.commutator_basis.len(), 2];
        let bch_commutator_basis_stride = compact_strides(&bch_commutator_basis_shape);

        let mut commutation_index_map = HashMap::new();

        for (i, commutator_term) in log_signature.bch_series.commutator_basis.iter().enumerate() {
            match commutator_term {
                t @ CommutatorTerm::Atom { atom, .. } => {
                    bch_commutator_basis.push([i32::from(*atom), -1]);
                    commutation_index_map.insert(t.atom_hash(), i);
                }
                t @ CommutatorTerm::Expression { left, right, .. } => {
                    let left_index = commutation_index_map[&left.atom_hash()] as i32;
                    let right_index = commutation_index_map[&right.atom_hash()] as i32;
                    bch_commutator_basis.push([left_index, right_index]);
                    commutation_index_map.insert(t.atom_hash(), i);
                }
            }
        }

        let bch_commutator_basis_flat = bch_commutator_basis
            .into_iter()
            .flatten()
            .collect::<Vec<_>>();

        let bch_commutator_basis_data = client.create(i32::as_bytes(&bch_commutator_basis_flat));

        Self {
            lie_series_data,
            bch_coefficients_data,
            bch_coefficients_len: log_signature.bch_series.coefficients.len(),
            bch_commutator_basis_data,
            bch_commutator_basis_shape,
            bch_commutator_basis_stride,
        }
    }

    pub fn empty(
        client: &ComputeClient<R::Server, R::Channel>,
        log_signature: &LogSignature<u8, T>,
    ) -> Self {
        let lie_series_data = LieSeriesCLData::<R, T>::empty(client, &log_signature.series);
        let bch_coefficients_data =
            client.create(T::as_bytes(&log_signature.bch_series.coefficients));
        let mut bch_commutator_basis =
            Vec::with_capacity(log_signature.bch_series.commutator_basis.len());

        let bch_commutator_basis_shape = vec![log_signature.bch_series.commutator_basis.len(), 2];
        let bch_commutator_basis_stride = compact_strides(&bch_commutator_basis_shape);

        let mut commutation_index_map = HashMap::new();

        for (i, commutator_term) in log_signature.bch_series.commutator_basis.iter().enumerate() {
            match commutator_term {
                t @ CommutatorTerm::Atom { atom, .. } => {
                    bch_commutator_basis.push([i32::from(*atom), -1]);
                    commutation_index_map.insert(t.atom_hash(), i);
                }
                t @ CommutatorTerm::Expression { left, right, .. } => {
                    let left_index = commutation_index_map[&left.atom_hash()] as i32;
                    let right_index = commutation_index_map[&right.atom_hash()] as i32;
                    bch_commutator_basis.push([left_index, right_index]);
                    commutation_index_map.insert(t.atom_hash(), i);
                }
            }
        }

        let bch_commutator_basis_flat = bch_commutator_basis
            .into_iter()
            .flatten()
            .collect::<Vec<_>>();

        let bch_commutator_basis_data = client.create(i32::as_bytes(&bch_commutator_basis_flat));

        Self {
            lie_series_data,
            bch_coefficients_data,
            bch_coefficients_len: log_signature.bch_series.coefficients.len(),
            bch_commutator_basis_data,
            bch_commutator_basis_shape,
            bch_commutator_basis_stride,
        }
    }

    #[must_use]
    pub fn kernel_arg<F: CubePrimitive>(&self) -> LogSignatureCLLaunch<'_, F, R> {
        let lie_series = self.lie_series_data.kernel_arg();
        let bch_coefficients = unsafe {
            ArrayArg::from_raw_parts::<T>(&self.bch_coefficients_data, self.bch_coefficients_len, 1)
        };
        let bch_commutator_basis = unsafe {
            TensorArg::from_raw_parts::<u32>(
                &self.bch_commutator_basis_data,
                &self.bch_commutator_basis_stride,
                &self.bch_commutator_basis_shape,
                1,
            )
        };
        LogSignatureCLLaunch::new(lie_series, bch_coefficients, bch_commutator_basis)
    }

    pub fn read_coefficients(self, client: &ComputeClient<R::Server, R::Channel>) -> Vec<T> {
        self.lie_series_data.read_coefficients(client)
    }
}

#[derive(CubeType, CubeLaunch)]
pub struct LogSignatureCL<T: CubePrimitive> {
    pub lie_series: LieSeriesCL<T>,
    pub bch_coefficients: Array<T>,
    pub bch_commutator_basis: Tensor<u32>,
}
