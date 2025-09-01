use std::{cmp::min, collections::HashMap};

use cubecl::{future::block_on, prelude::*};
use lie_cl::lie_series::{
    LieSeriesCLData, lie_series_addition, lie_series_commutation, lie_series_scalar_multiplication,
};
use lyndon_rs::{LyndonBasis, Sort};
use signature_rs::{CommutatorTerm, LogSignature};

const MAX_THREADS: u32 = 1024;
const LINE_SIZE: u8 = 4;

pub struct LogSignatureCLData<'a, R, T, U> {
    /// The log signature structure
    log_signature: &'a LogSignature<u8, U>,
    /// The necessary information for calculating Lie Series commutations
    lie_series_data: LieSeriesCLData<R, T>,
    /// The dependencies for each term in the series, e.g.
    /// `A = [0, -1]`, `B = [1, -1]`, `[A, B] = [0, 1]`, `[A, [A, B]] = [0, 2]`, ...
    bch_commutator_basis: Vec<[i32; 2]>,
}

impl<R: Runtime, T: Numeric + CubeElement + From<U>, U> Clone for LogSignatureCLData<'_, R, T, U> {
    fn clone(&self) -> Self {
        Self {
            lie_series_data: self.lie_series_data.clone(),
            log_signature: self.log_signature,
            bch_commutator_basis: self.bch_commutator_basis.clone(),
        }
    }
}

impl<'a, R: Runtime, T: Copy + Numeric + CubeElement + From<U>, U: Copy>
    LogSignatureCLData<'a, R, T, U>
{
    #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
    pub fn new(
        client: &ComputeClient<R::Server, R::Channel>,
        log_signature: &'a LogSignature<u8, U>,
    ) -> Self {
        let lie_series_data = LieSeriesCLData::<R, T>::new(client, &log_signature.series);

        let mut bch_commutator_basis =
            Vec::with_capacity(log_signature.bch_series.commutator_basis.len());

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

        Self {
            log_signature,
            lie_series_data,
            bch_commutator_basis,
        }
    }

    #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
    pub fn empty(
        client: &ComputeClient<R::Server, R::Channel>,
        log_signature: &'a LogSignature<u8, U>,
    ) -> Self {
        let lie_series_data = LieSeriesCLData::<R, T>::empty(client, &log_signature.series);
        let mut bch_commutator_basis =
            Vec::with_capacity(log_signature.bch_series.commutator_basis.len());

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

        Self {
            log_signature,
            lie_series_data,
            bch_commutator_basis,
        }
    }

    #[must_use]
    #[allow(
        clippy::too_many_lines,
        clippy::cast_possible_wrap,
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::similar_names
    )]
    pub fn concat(&self, client: &ComputeClient<R::Server, R::Channel>, other: &Self) -> Self {
        let basis = LyndonBasis::<u8>::new(2, Sort::Lexicographical);
        let terms_per_degree =
            basis.number_of_words_per_degree(self.log_signature.bch_series.max_degree);
        let mut commutation_terms = Vec::with_capacity(self.bch_commutator_basis.len());
        commutation_terms.push(self.lie_series_data.clone());
        commutation_terms.push(other.lie_series_data.clone());
        for _ in 2..self.bch_commutator_basis.len() {
            commutation_terms.push(LieSeriesCLData::<R, T>::empty(
                client,
                &self.log_signature.series,
            ));
        }

        let dim_i = self.lie_series_data.basis_map_shape[0] as u32;
        let dim_j = self.lie_series_data.basis_map_shape[1] as u32;
        let dim_k = self.lie_series_data.basis_map_shape[2] as u32;

        let cube_x = min(dim_i, 32);
        let cube_y = min(dim_j, 32);
        let cube_z = min(dim_k, 1024 / (cube_x * cube_y));

        let cubes_x = dim_i.div_ceil(cube_x);
        let cubes_y = dim_j.div_ceil(cube_y);
        let cubes_z = dim_k.div_ceil(cube_z);

        let cube_count = CubeCount::Static(cubes_x, cubes_y, cubes_z);
        let cube_dim = CubeDim::new(cube_x, cube_y, cube_z);

        let mut cumsum = 0;
        for num_terms in terms_per_degree {
            for i in 0..num_terms {
                let plan = self.bch_commutator_basis[i + cumsum];
                if plan[1] == -1 {
                    continue;
                }
                let left_index = plan[0] as usize;
                let right_index = plan[1] as usize;

                let left_term = &commutation_terms[left_index];
                let right_term = &commutation_terms[right_index];
                let output_term = &commutation_terms[i + cumsum];

                unsafe {
                    lie_series_commutation::launch_unchecked::<T, R>(
                        client,
                        cube_count.clone(),
                        cube_dim,
                        left_term.kernel_arg(1),
                        right_term.kernel_arg(1),
                        output_term.kernel_arg(1),
                    );
                }
            }
            cumsum += num_terms;
            // Ensure the previous level's calculations are complete.
            block_on(client.sync());
        }

        let num_blocks = self
            .log_signature
            .bch_series
            .coefficients
            .len()
            .div_ceil(MAX_THREADS as usize) as u32;

        let cube_count = CubeCount::Static(num_blocks, 1, 1);
        let cube_dim = CubeDim::new(MAX_THREADS, 1, 1);

        for (coefficient, term) in self
            .log_signature
            .bch_series
            .coefficients
            .iter()
            .zip(commutation_terms.iter())
        {
            unsafe {
                lie_series_scalar_multiplication::launch_unchecked::<T, R>(
                    client,
                    cube_count.clone(),
                    cube_dim,
                    term.kernel_arg(LINE_SIZE),
                    ScalarArg::new(<T as std::convert::From<U>>::from(*coefficient)),
                    term.kernel_arg(LINE_SIZE),
                );
            }
        }
        block_on(client.sync());

        let reduction = &commutation_terms[0];

        for i in self
            .log_signature
            .bch_series
            .coefficients
            .iter()
            .enumerate()
            .filter(|(_, c)| <T as std::convert::From<U>>::from(**c) != T::from_int(0))
            .map(|(x, _)| x)
            .skip(1)
        {
            unsafe {
                lie_series_addition::launch_unchecked::<T, R>(
                    client,
                    cube_count.clone(),
                    cube_dim,
                    reduction.kernel_arg(LINE_SIZE),
                    commutation_terms[i].kernel_arg(LINE_SIZE),
                    reduction.kernel_arg(LINE_SIZE),
                );
            }
            block_on(client.sync());
        }

        Self {
            log_signature: self.log_signature,
            lie_series_data: reduction.clone(),
            bch_commutator_basis: self.bch_commutator_basis.clone(),
        }
    }

    pub fn read_coefficients(self, client: &ComputeClient<R::Server, R::Channel>) -> Vec<T> {
        self.lie_series_data.read_coefficients(client)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use cubecl::cuda::{CudaDevice, CudaRuntime};
    use ordered_float::NotNan;
    use rstest::rstest;
    use signature_rs::LogSignatureBuilder;

    #[rstest]
    fn test_log_sig_concat() {
        let builder = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(2)
            .with_max_degree(3);
        let mut a = builder.build();
        let mut b = builder.build();
        a.series.coefficients = [1., 2., 3., 4., 5.]
            .map(|x| unsafe { NotNan::new_unchecked(x) })
            .to_vec();
        b.series.coefficients = [6., 7., 8., 9., 10.]
            .map(|x| unsafe { NotNan::new_unchecked(x) })
            .to_vec();

        let expected_concat = a
            .concatenate(&b)
            .series
            .coefficients
            .into_iter()
            .map(ordered_float::NotNan::into_inner)
            .collect::<Vec<_>>();

        let client = CudaRuntime::client(&CudaDevice::default());
        let log_sig_a = LogSignatureCLData::<CudaRuntime, f32, NotNan<f32>>::new(&client, &a);
        let log_sig_b = LogSignatureCLData::<CudaRuntime, f32, NotNan<f32>>::new(&client, &b);

        let concat = log_sig_a.concat(&client, &log_sig_b);
        let coefficients = concat.read_coefficients(&client);

        dbg!(&coefficients);
        for (c, e_c) in coefficients.into_iter().zip(expected_concat) {
            assert!((c - e_c).abs() < 0.001, "{c} != {e_c}");
        }
    }
}
