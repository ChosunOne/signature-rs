//! Basic log signature computation example.
//!
//! This example demonstrates how to compute log signatures from 2D path data
//! using the signature-rs library.

use ndarray::array;
use num_traits::{Signed, ToPrimitive};
use ordered_float::NotNan;
use signature_rs::prelude::*;

fn main() {
    println!("=== Log Signature Computation Example ===\n");

    // Create a simple 2D path: a triangle
    let path = array![
        [0.0, 0.0], // Start at origin
        [1.0, 0.0], // Move right
        [0.5, 1.0], // Move up and left
        [0.0, 0.0]  // Return to origin
    ];

    println!("Input path (triangle):");
    for (i, point) in path.rows().into_iter().enumerate() {
        println!("  Point {}: [{:.1}, {:.1}]", i, point[0], point[1]);
    }
    println!();

    // Convert to NotNan for numerical stability
    let path = path.mapv(|v| NotNan::new(v).expect("value to be a number"));

    // Create a log signature builder
    let builder = LogSignatureBuilder::<ENotation>::new()
        .with_num_dimensions(2)
        .with_max_degree(3);

    println!("Builder configuration:");
    println!("  Dimensions: {}", 2);
    println!("  Max degree: {}", 3);
    println!();

    // Compute the log signature
    let log_sig = builder.build_from_path(&path.view());

    println!("Log signature computation results:");
    println!("  Number of basis terms: {}", log_sig.series.basis.len());
    println!(
        "  Number of coefficients: {}",
        log_sig.series.coefficients.len()
    );
    println!();

    // Display the first few terms
    println!("First few log signature terms:");
    for (i, (basis_term, coeff)) in log_sig
        .series
        .basis
        .iter()
        .zip(log_sig.series.coefficients.iter())
        .take(10)
        .enumerate()
    {
        println!("  Term {}: {} = {:.6}", i + 1, basis_term, coeff);
    }

    if log_sig.series.basis.len() > 10 {
        println!("  ... and {} more terms", log_sig.series.basis.len() - 10);
    }
    println!();

    // Demonstrate path signature properties
    println!("=== Path Signature Properties ===");

    // Compute signatures for path segments
    let segment1 = path.slice(ndarray::s![0..=2, ..]);
    let segment2 = path.slice(ndarray::s![2.., ..]);

    let log_sig1 = builder.build_from_path(&segment1);
    let log_sig2 = builder.build_from_path(&segment2);

    println!(
        "Segment 1 (first 3 points): {} coefficients",
        log_sig1.series.coefficients.len()
    );
    println!(
        "Segment 2 (last 2 points): {} coefficients",
        log_sig2.series.coefficients.len()
    );

    // Concatenate the segments
    let concatenated = log_sig1.concatenate(&log_sig2);
    println!(
        "Concatenated signature: {} coefficients",
        concatenated.series.coefficients.len()
    );

    // Compare with full path signature
    let diff: f64 = log_sig
        .series
        .coefficients
        .iter()
        .zip(concatenated.series.coefficients.iter())
        .map(|(a, b)| {
            let diff = *a - *b;
            diff.abs().to_f64().unwrap_or(0.0)
        })
        .sum();

    println!("Difference from full path: {:.2e}", diff);

    if diff < 1e-10 {
        println!("✓ Concatenation property verified!");
    } else {
        println!("⚠ Unexpected difference - this may indicate numerical issues");
    }

    println!("\n=== Example Complete ===");
}
