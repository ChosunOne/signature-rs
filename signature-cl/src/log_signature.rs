use std::collections::HashMap;

use cubecl::{prelude::*, server::Handle};
use lie_cl::lie_series::LieSeriesCLData;
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
    bch_commutator_basis: Vec<[i32; 2]>,
}

impl<R: Runtime, T: Numeric + CubeElement> Clone for LogSignatureCLData<R, T> {
    fn clone(&self) -> Self {
        Self {
            lie_series_data: self.lie_series_data.clone(),
            bch_coefficients_data: self.bch_coefficients_data.clone(),
            bch_coefficients_len: self.bch_coefficients_len,
            bch_commutator_basis: self.bch_commutator_basis.clone(),
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
            lie_series_data,
            bch_coefficients_data,
            bch_coefficients_len: log_signature.bch_series.coefficients.len(),
            bch_commutator_basis,
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
            lie_series_data,
            bch_coefficients_data,
            bch_coefficients_len: log_signature.bch_series.coefficients.len(),
            bch_commutator_basis,
        }
    }

    pub fn concat(&self, other: &Self) {
        todo!()
    }

    pub fn read_coefficients(self, client: &ComputeClient<R::Server, R::Channel>) -> Vec<T> {
        self.lie_series_data.read_coefficients(client)
    }
}

#[cfg(test)]
mod test {
    use std::collections::HashMap;

    use lie_rs::{BchSeriesGenerator, LieSeries, LieSeriesGenerator};
    use lyndon_rs::{LyndonBasis, Sort};
    use ordered_float::NotNan;
    use signature_rs::CommutatorTerm;

    #[test]
    fn test_bch_series_plan() {
        let basis = LyndonBasis::<u8>::new(2, Sort::Lexicographical);
        let bch_series: LieSeries<u8, NotNan<f32>> =
            BchSeriesGenerator::new(basis, 5).generate_lie_series();

        let mut commutation_index_map = HashMap::new();
        let mut bch_commutator_basis = Vec::with_capacity(bch_series.commutator_basis.len());

        for (i, commutator_term) in bch_series.commutator_basis.iter().enumerate() {
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

        println!("Index \t Commutator \t Degree");
        for (i, plan) in bch_commutator_basis.iter().enumerate() {
            println!(
                "[{i}]:\t [{}, {}] \t{}",
                plan[0],
                plan[1],
                &bch_series.commutator_basis[i].degree()
            );
        }
        todo!()
    }
}
