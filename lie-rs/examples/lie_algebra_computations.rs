//! Free Lie algebra computations and properties example.
//!
//! This example explores the structure of free Lie algebras using
//! the tools provided by lie-rs and lyndon-rs.

use lie_rs::prelude::*;
use lyndon_rs::prelude::*;
use num_rational::Ratio;
use num_traits::{Signed, Zero};
use std::collections::HashMap;

fn main() {
    println!("=== Free Lie Algebra Computations ===\n");

    // Explore different alphabet sizes
    for alphabet_size in 2..=4 {
        println!("=== Free Lie Algebra on {} Generators ===", alphabet_size);
        explore_free_lie_algebra(alphabet_size);
        println!();
    }

    println!("=== Bernoulli Numbers and BCH ===");
    demonstrate_bernoulli_connection();

    println!("\n=== Lie Series Operations ===");
    demonstrate_lie_series_operations();

    println!("\n=== Example Complete ===");
}

fn explore_free_lie_algebra(alphabet_size: usize) {
    let basis = LyndonBasis::<ENotation>::new(alphabet_size, Sort::Lexicographical);
    let max_degree = 6;

    // Count dimensions by degree
    let counts = basis.number_of_words_per_degree(max_degree);
    println!("Dimensions of graded components:");

    let mut total_dim = 0;
    for (degree, count) in counts.iter().enumerate().skip(1) {
        if *count > 0 {
            total_dim += count;
            println!("  deg {}: {} (total: {})", degree, count, total_dim);
        }
    }

    // Generate basis elements up to small degree for display
    let words = basis.generate_basis(4);
    println!("\nBasis elements up to degree 4:");

    let mut words_by_degree = HashMap::new();
    for word in &words {
        words_by_degree
            .entry(word.letters.len())
            .or_insert(Vec::new())
            .push(word);
    }

    for degree in 1..=4 {
        if let Some(degree_words) = words_by_degree.get(&degree) {
            println!("  Degree {}: {}", degree, degree_words.len());
            for (i, word) in degree_words.iter().enumerate() {
                if i < 8 {
                    // Show first 8 elements
                    println!("    {}", word);
                }
            }
            if degree_words.len() > 8 {
                println!("    ... and {} more", degree_words.len() - 8);
            }
        }
    }

    // Show growth rate
    let growth_info = analyze_growth_rate(&counts, alphabet_size);
    println!("\nGrowth analysis:");
    println!("  {}", growth_info);
}

fn analyze_growth_rate(counts: &[usize], alphabet_size: usize) -> String {
    // The dimension grows roughly like n^(alphabet_size-1) / n! for large n
    // For 2 generators, it's roughly exponential with base related to alphabet size

    if counts.len() < 4 {
        return "Insufficient data for growth analysis".to_string();
    }

    let non_zero_counts: Vec<_> = counts
        .iter()
        .enumerate()
        .skip(1)
        .filter(|&(_, &count)| count > 0)
        .map(|(deg, &count)| (deg, count))
        .collect();

    if non_zero_counts.len() < 3 {
        return "Need more degrees for analysis".to_string();
    }

    // Look at ratios between consecutive terms
    let mut ratios = Vec::new();
    for i in 1..non_zero_counts.len() {
        let curr = non_zero_counts[i].1 as f64;
        let prev = non_zero_counts[i - 1].1 as f64;
        if prev > 0.0 {
            ratios.push(curr / prev);
        }
    }

    if ratios.is_empty() {
        return "Cannot compute growth ratios".to_string();
    }

    let avg_ratio = ratios.iter().sum::<f64>() / ratios.len() as f64;

    format!(
        "Average growth ratio: {:.2} (base ≈ {}^(1/n))",
        avg_ratio, alphabet_size
    )
}

fn demonstrate_bernoulli_connection() {
    println!("Connection between Bernoulli numbers and BCH coefficients:");
    println!();

    // Generate some Bernoulli numbers
    // Generate first few Bernoulli numbers manually for demonstration
    let bernoulli_nums = vec![
        Ratio::new(1, 1),   // B₀ = 1
        Ratio::new(-1, 2),  // B₁ = -1/2
        Ratio::new(1, 6),   // B₂ = 1/6
        Ratio::new(0, 1),   // B₃ = 0
        Ratio::new(-1, 30), // B₄ = -1/30
        Ratio::new(0, 1),   // B₅ = 0
        Ratio::new(1, 42),  // B₆ = 1/42
        Ratio::new(0, 1),   // B₇ = 0
        Ratio::new(-1, 30), // B₈ = -1/30
        Ratio::new(0, 1),   // B₉ = 0
    ];

    println!("First few Bernoulli numbers B_n:");
    for (n, b_n) in bernoulli_nums.iter().enumerate().take(8) {
        if n == 1 {
            println!("  B_{} = {} (using positive convention)", n, b_n);
        } else {
            println!("  B_{} = {}", n, b_n);
        }
    }

    println!();
    println!("In the BCH series, Bernoulli numbers appear in the coefficients");
    println!("of certain higher-order terms. For example:");
    println!("- The coefficient of [X,[X,[X,Y]]] involves B_2 = 1/6");
    println!("- Higher order terms involve higher Bernoulli numbers");
    println!();

    // Show the connection more explicitly
    println!("BCH coefficient structure:");
    println!("- Degree 2: involves 1/2 (related to B_1)");
    println!("- Degree 3: involves 1/12 (related to factorials)");
    println!("- Degree 4: involves 1/24, -1/12 (related to B_2 = 1/6)");
    println!("- Higher degrees: more complex Bernoulli combinations");
}

fn demonstrate_lie_series_operations() {
    println!("Demonstrating Lie series operations:");
    println!();

    // Create two small BCH series
    let basis = LyndonBasis::<ENotation>::new(2, Sort::Lexicographical);
    let generator1 = BchSeriesGenerator::<ENotation>::new(basis.clone(), 3);
    let generator2 = BchSeriesGenerator::<ENotation>::new(basis.clone(), 3);

    let series1: LieSeries<ENotation, Ratio<i64>> = generator1.generate_lie_series();
    let series2: LieSeries<ENotation, Ratio<i64>> = generator2.generate_lie_series();

    println!("Series 1 (BCH up to degree 3):");
    display_series_summary(&series1);

    println!("\nSeries 2 (identical series for demonstration):");
    display_series_summary(&series2);

    // Show series properties
    println!("\nSeries properties:");
    analyze_series_properties(&series1);

    // Demonstrate coefficient access
    println!("\nCoefficient extraction:");
    demonstrate_coefficient_access(&series1);
}

fn display_series_summary(series: &LieSeries<ENotation, Ratio<i64>>) {
    let non_zero_terms = series
        .coefficients
        .iter()
        .filter(|&coeff| !coeff.is_zero())
        .count();

    println!("  Total terms: {}", series.coefficients.len());
    println!("  Non-zero terms: {}", non_zero_terms);

    // Show first few non-zero terms
    let mut shown = 0;
    for (basis_elem, coeff) in series.basis.iter().zip(series.coefficients.iter()) {
        if !coeff.is_zero() && shown < 5 {
            println!("    {} * {}", format_rational(coeff), basis_elem);
            shown += 1;
        }
    }
    if non_zero_terms > 5 {
        println!("    ... and {} more terms", non_zero_terms - 5);
    }
}

fn analyze_series_properties(series: &LieSeries<ENotation, Ratio<i64>>) {
    // Count terms by degree
    let mut degree_counts = HashMap::new();
    for (basis_elem, coeff) in series.basis.iter().zip(series.coefficients.iter()) {
        if !coeff.is_zero() {
            let degree = basis_elem.letters.len();
            *degree_counts.entry(degree).or_insert(0) += 1;
        }
    }

    println!("  Terms by degree:");
    for degree in 1..=6 {
        if let Some(count) = degree_counts.get(&degree) {
            println!("    Degree {}: {} terms", degree, count);
        }
    }

    // Compute total "magnitude" (sum of absolute coefficients)
    let total_magnitude: Ratio<i64> = series.coefficients.iter().map(|c| c.abs()).sum();

    println!("  Total magnitude: {}", format_rational(&total_magnitude));
}

fn demonstrate_coefficient_access(series: &LieSeries<ENotation, Ratio<i64>>) {
    println!("Accessing specific coefficients:");

    // Look for linear terms (degree 1)
    for (basis_elem, coeff) in series.basis.iter().zip(series.coefficients.iter()) {
        if basis_elem.letters.len() == 1 && !coeff.is_zero() {
            println!(
                "  Linear term {}: coefficient = {}",
                basis_elem,
                format_rational(coeff)
            );
        }
    }

    // Look for degree 2 terms
    for (basis_elem, coeff) in series.basis.iter().zip(series.coefficients.iter()) {
        if basis_elem.letters.len() == 2 && !coeff.is_zero() {
            println!(
                "  Quadratic term {}: coefficient = {}",
                basis_elem,
                format_rational(coeff)
            );
        }
    }
}

fn format_rational(r: &Ratio<i64>) -> String {
    if *r.denom() == 1 {
        format!("{}", r.numer())
    } else {
        format!("{}/{}", r.numer(), r.denom())
    }
}
