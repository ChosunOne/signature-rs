//! Baker-Campbell-Hausdorff series computation example.
//!
//! This example demonstrates how to compute the Baker-Campbell-Hausdorff series
//! using the lie-rs library, showing the expansion of log(exp(X)exp(Y)).

use lie_rs::prelude::*;
use lyndon_rs::prelude::*;
use num_rational::Ratio;
use num_traits::{Signed, Zero};

fn main() {
    println!("=== Baker-Campbell-Hausdorff Series Example ===\n");

    // Create a Lyndon basis for 2 generators (X and Y)
    let basis = LyndonBasis::<ENotation>::new(2, Sort::Lexicographical);
    println!(
        "Created Lyndon basis with {} generators",
        basis.alphabet_size
    );

    // Show the computation for different degrees
    for max_degree in 2..=5 {
        println!("\n=== BCH Series up to degree {} ===", max_degree);
        demonstrate_bch_series(&basis, max_degree);
    }

    println!("\n=== Detailed Analysis: Degree 5 ===");
    detailed_analysis(&basis, 5);

    println!("\n=== BCH Series Properties ===");
    explore_bch_properties();

    println!("\n=== Example Complete ===");
}

fn demonstrate_bch_series(basis: &LyndonBasis<ENotation>, max_degree: usize) {
    // Create BCH series generator
    let generator = BchSeriesGenerator::<ENotation>::new(basis.clone(), max_degree);

    // Generate the Lie series
    let bch_series: LieSeries<ENotation, Ratio<i64>> = generator.generate_lie_series();

    println!("BCH series: log(exp(X)exp(Y)) = X + Y + [X,Y]/2 + ...");
    println!(
        "Generated series with {} terms",
        bch_series.coefficients.len()
    );

    // Display the series terms
    println!("\nSeries expansion:");
    for (i, (basis_element, coeff)) in bch_series
        .basis
        .iter()
        .zip(bch_series.coefficients.iter())
        .enumerate()
    {
        if !coeff.is_zero() {
            let coeff_str = format_coefficient(coeff);
            println!("  Term {}: {} * {}", i + 1, coeff_str, basis_element);
        }
    }

    // Show series grouped by degree
    println!("\nTerms grouped by degree:");
    let mut terms_by_degree: std::collections::HashMap<usize, Vec<(String, String)>> =
        std::collections::HashMap::new();

    for (basis_element, coeff) in bch_series.basis.iter().zip(bch_series.coefficients.iter()) {
        if !coeff.is_zero() {
            let degree = basis_element.letters.len();
            let coeff_str = format_coefficient(coeff);
            terms_by_degree
                .entry(degree)
                .or_default()
                .push((coeff_str, format!("{}", basis_element)));
        }
    }

    for degree in 1..=max_degree {
        if let Some(terms) = terms_by_degree.get(&degree) {
            println!("  Degree {}: {} terms", degree, terms.len());
            for (coeff, basis_elem) in terms {
                println!("    {} * {}", coeff, basis_elem);
            }
        }
    }
}

fn detailed_analysis(basis: &LyndonBasis<ENotation>, max_degree: usize) {
    let generator = BchSeriesGenerator::<ENotation>::new(basis.clone(), max_degree);
    let bch_series: LieSeries<ENotation, Ratio<i64>> = generator.generate_lie_series();

    println!("Detailed analysis of BCH series coefficients:");
    println!();

    // Analyze coefficient patterns
    let mut degree_stats = std::collections::HashMap::new();

    for (basis_element, coeff) in bch_series.basis.iter().zip(bch_series.coefficients.iter()) {
        if !coeff.is_zero() {
            let degree = basis_element.letters.len();
            let entry = degree_stats.entry(degree).or_insert((0, Vec::new()));
            entry.0 += 1;
            entry.1.push(coeff.clone());
        }
    }

    for degree in 1..=max_degree {
        if let Some((count, coeffs)) = degree_stats.get(&degree) {
            println!("Degree {} analysis:", degree);
            println!("  Number of terms: {}", count);

            // Find common denominators
            let denominators: Vec<i64> = coeffs.iter().map(|c| *c.denom()).collect();
            let max_denom = denominators.iter().max().unwrap_or(&1);
            println!("  Largest denominator: {}", max_denom);

            // Show coefficient magnitudes
            let max_coeff = coeffs
                .iter()
                .map(|c| c.abs())
                .max()
                .unwrap_or(Ratio::new(0, 1));
            println!(
                "  Largest coefficient magnitude: {}",
                format_coefficient(&max_coeff)
            );
            println!();
        }
    }

    // Show the classical BCH formula components
    println!("Classical BCH formula verification:");
    println!("log(exp(X)exp(Y)) = X + Y + (1/2)[X,Y] + (1/12)([X,[X,Y]] + [Y,[Y,X]]) + ...");
    println!();

    // Verify known coefficients
    verify_known_coefficients(&bch_series);
}

fn explore_bch_properties() {
    println!("Properties of the Baker-Campbell-Hausdorff series:");
    println!();

    println!("1. Convergence: The BCH series converges when ||X|| + ||Y|| < log(2)");
    println!("   for matrices in certain norms.");
    println!();

    println!("2. Symmetry: BCH(X,Y) = -BCH(-Y,-X)");
    println!("   This reflects the group property.");
    println!();

    println!("3. Degree structure: Terms of degree n appear with denominators");
    println!("   that are divisors of n! (but typically much smaller).");
    println!();

    println!("4. Lyndon basis: Each term corresponds to a Lyndon word,");
    println!("   providing a natural basis for the free Lie algebra.");
    println!();

    println!("5. Applications:");
    println!("   - Solving differential equations");
    println!("   - Magnus expansion in physics");
    println!("   - Geometric integration methods");
    println!("   - Quantum mechanics (time evolution)");
}

fn format_coefficient(coeff: &Ratio<i64>) -> String {
    if *coeff.denom() == 1 {
        format!("{}", coeff.numer())
    } else {
        format!("{}/{}", coeff.numer(), coeff.denom())
    }
}

fn verify_known_coefficients(bch_series: &LieSeries<ENotation, Ratio<i64>>) {
    println!("Verifying known BCH coefficients:");

    // Look for specific terms and check their coefficients
    for (i, (basis_element, coeff)) in bch_series
        .basis
        .iter()
        .zip(bch_series.coefficients.iter())
        .enumerate()
    {
        if !coeff.is_zero() {
            let generators = &basis_element.letters;

            // Check for specific patterns
            match generators.len() {
                1 => {
                    // X or Y terms should have coefficient 1
                    if *coeff == Ratio::new(1, 1) {
                        println!("  ✓ Linear term {} has coefficient 1", basis_element);
                    }
                }
                2 => {
                    // [X,Y] should have coefficient 1/2
                    // Check for [X,Y] pattern - this would be e₁e₂
                    if generators.len() == 2 {
                        let expected = Ratio::new(1, 2);
                        if *coeff == expected {
                            println!("  ✓ [X,Y] has coefficient 1/2");
                        } else {
                            println!(
                                "  ✗ [X,Y] has coefficient {} (expected 1/2)",
                                format_coefficient(coeff)
                            );
                        }
                    }
                }
                3 => {
                    // Check for known degree 3 terms
                    if *coeff == Ratio::new(1, 12) || *coeff == Ratio::new(-1, 12) {
                        println!(
                            "  ✓ Degree 3 term {} has coefficient {}",
                            basis_element,
                            format_coefficient(coeff)
                        );
                    }
                }
                _ => {
                    // Just show a few higher degree terms
                    if i < 20 {
                        println!(
                            "  → {} has coefficient {}",
                            basis_element,
                            format_coefficient(coeff)
                        );
                    }
                }
            }
        }
    }
}
