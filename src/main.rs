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

    println!("Non-basis terms");
    for term in &non_basis_terms {
        println!("{term}");
    }
    println!();

    let mut non_basis_terms_queue = non_basis_terms.clone();
    let mut normalized_terms = vec![];
    let mut processed_terms = HashSet::new();
    while let Some(term) = non_basis_terms_queue.pop() {
        if processed_terms.contains(&term) {
            println!("Skipping {term}");
            println!();
            continue;
        }
        processed_terms.insert(term.clone());
        println!("{term}");

        let Some((left_term, right_term)) = term.jacobi_identity() else {
            continue;
        };
        println!("lj: {left_term}");
        println!("rj: {right_term}");

        let left_term_key = match &left_term {
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

        let mut left_normalized = false;
        let mut right_normalized = false;

        if series_a.commutator_basis_map.contains_key(&left_term_key) {
            normalized_terms.push(left_term);
            left_normalized = true;
        } else if !left_term.is_zero() {
            non_basis_terms_queue.push(left_term);
        } else {
            left_normalized = true;
        }

        let right_term_key = match &right_term {
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

        if series_a.commutator_basis_map.contains_key(&right_term_key) {
            normalized_terms.push(right_term);
            right_normalized = true;
        } else if !right_term.is_zero() {
            non_basis_terms_queue.push(right_term);
        } else {
            right_normalized = true;
        }

        if left_normalized && right_normalized {
            processed_terms.remove(&term);
        }
        println!("left_normalized: {left_normalized}");
        println!("right_normalized: {right_normalized}");
        println!();
    }

    println!("Normalized Non-Basis Terms");
    for term in &normalized_terms {
        println!("{term}");
        processed_terms.remove(term);
    }
    println!();

    for term in &non_basis_terms {
        processed_terms.remove(term);
    }

    println!("Terms that failed to normalize");
    for term in processed_terms {
        match term {
            CommutatorTerm::Atom { coefficient, .. }
            | CommutatorTerm::Expression { coefficient, .. } => {
                if coefficient == 0 {
                    continue;
                }
            }
        }
        println!("{term}");
    }
    println!();
}
