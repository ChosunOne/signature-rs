use std::collections::HashSet;

use signature_rs::comm;
use signature_rs::commutator::{Commutator, CommutatorTerm, FormalIndeterminate};
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
            if !series_a.commutator_basis_index_map.contains_key(&term_key) {
                non_basis_terms.push(term.clone());
            }
        }
    }

    println!("Non-basis terms");
    for (i, term) in non_basis_terms.iter().enumerate() {
        println!("{i}: {term}");
    }
    println!();

    println!("Non-basis Commutator Basis Map");
    for (i, non_basis_term) in non_basis_terms.iter().enumerate() {
        let non_basis_term_key = match non_basis_term {
            CommutatorTerm::Atom { atom, .. } => CommutatorTerm::Atom {
                coefficient: 1,
                atom: *atom,
            },
            CommutatorTerm::Expression { left, right, .. } => CommutatorTerm::Expression {
                coefficient: 1,
                left: left.clone(),
                right: right.clone(),
            },
        };

        let v = &series_a.commutator_basis_map[&non_basis_term_key];

        print!("{i}: [");
        for basis_term in v {
            print!("{basis_term}, ");
        }
        print!("]");
        println!();

        // println!("Indeterminate expansion");
        // for basis_term in v {
        //     println!("{basis_term}");
        //     let indeterminates = Vec::<FormalIndeterminate<ENotation, i32>>::from(basis_term);
        //     for indeterminate in &indeterminates {
        //         println!("{indeterminate}");
        //     }
        //     println!();
        // }
        //
        // println!("Non-basis indeterminate expansion");
        // let indeterminates = Vec::<FormalIndeterminate<ENotation, i32>>::from(&non_basis_term_key);
        // for indeterminate in &indeterminates {
        //     println!("{indeterminate}");
        // }
        // println!();
    }
}
