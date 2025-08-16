use std::collections::HashSet;

use signature_rs::comm;
use signature_rs::commutator::{Commutator, CommutatorTerm};
use signature_rs::generators::ENotation;
use signature_rs::lie_series::LieSeries;
use signature_rs::lyndon::{LyndonBasis, Sort};

#[allow(clippy::too_many_lines)]
pub fn main() {
    let lyndon_basis = LyndonBasis::<ENotation>::new(3, Sort::Lexicographical);
    let basis = lyndon_basis.generate_basis(4);
    let series_a = LieSeries::new(
        basis.clone(),
        vec![
            0, 0, 0, 0, 0, 0, -7, -14, 14, 14, 49, 42, -14, 21, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0,
        ],
    );
    let series_b = LieSeries::new(
        basis,
        vec![
            5, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0,
        ],
    );

    println!("Commutator Basis");
    for (i, cb) in series_a.commutator_basis.iter().enumerate() {
        println!("{i}: {cb}");
    }
    println!();

    let mut terms = vec![];
    for i in 0..series_a.basis.len() {
        for j in 0..series_b.basis.len() {
            if i == j {
                continue;
            }
            if series_a.coefficients[i] == 0 || series_b.coefficients[j] == 0 {
                continue;
            }
            let mut commutation = comm![series_a.commutator_basis[i], series_b.commutator_basis[j]]
                * series_a.coefficients[i]
                * series_b.coefficients[j];
            commutation.lyndon_sort();
            terms.push(commutation);
        }
    }

    terms.sort();
    println!("Multiplied terms");
    for (i, term) in terms.iter().enumerate() {
        println!("{i}: {term}");
    }
    println!();

    let mut non_basis_terms = vec![];
    for term in &terms {
        if let CommutatorTerm::Expression { left, right, .. } = term {
            let term_key = CommutatorTerm::Expression {
                coefficient: 1,
                left: left.clone(),
                right: right.clone(),
            };
            if !series_a.commutator_basis_map.contains_key(&term_key) {
                non_basis_terms.push(term.clone());
            }
        }
    }
}
