//! Time series analysis using log signatures.
//!
//! This example demonstrates how to use log signatures to analyze
//! and compare different time series patterns.

use num_traits::{Signed, ToPrimitive};
use ordered_float::NotNan;
use signature_rs::prelude::*;

fn main() {
    println!("=== Time Series Analysis with Log Signatures ===\n");

    // Create different time series patterns
    let linear_trend = create_linear_series();
    let oscillating = create_oscillating_series();
    let random_walk = create_random_walk();

    println!("Generated three time series:");
    println!("1. Linear trend: {} points", linear_trend.nrows());
    println!("2. Oscillating: {} points", oscillating.nrows());
    println!("3. Random walk: {} points", random_walk.nrows());
    println!();

    // Convert to NotNan
    let linear_trend = linear_trend.mapv(|v| NotNan::new(v).expect("value to be a number"));
    let oscillating = oscillating.mapv(|v| NotNan::new(v).expect("value to be a number"));
    let random_walk = random_walk.mapv(|v| NotNan::new(v).expect("value to be a number"));

    // Create signature builder
    let builder = LogSignatureBuilder::<ENotation>::new()
        .with_num_dimensions(2)
        .with_max_degree(4);

    // Compute log signatures
    let sig_linear = builder.build_from_path(&linear_trend.view());
    let sig_oscillating = builder.build_from_path(&oscillating.view());
    let sig_random = builder.build_from_path(&random_walk.view());

    println!("Log signature dimensions:");
    println!(
        "  Linear trend: {} terms",
        sig_linear.series.coefficients.len()
    );
    println!(
        "  Oscillating: {} terms",
        sig_oscillating.series.coefficients.len()
    );
    println!(
        "  Random walk: {} terms",
        sig_random.series.coefficients.len()
    );
    println!();

    // Analyze signature characteristics
    println!("=== Signature Analysis ===");

    analyze_signature("Linear Trend", &sig_linear);
    analyze_signature("Oscillating", &sig_oscillating);
    analyze_signature("Random Walk", &sig_random);

    // Compare signatures using L2 distance
    println!("=== Signature Comparisons (L2 Distance) ===");

    let dist_linear_osc = compute_l2_distance(&sig_linear, &sig_oscillating);
    let dist_linear_random = compute_l2_distance(&sig_linear, &sig_random);
    let dist_osc_random = compute_l2_distance(&sig_oscillating, &sig_random);

    println!("Linear vs Oscillating: {:.6}", dist_linear_osc);
    println!("Linear vs Random: {:.6}", dist_linear_random);
    println!("Oscillating vs Random: {:.6}", dist_osc_random);
    println!();

    // Find most similar pair
    let distances = [
        ("Linear vs Oscillating", dist_linear_osc),
        ("Linear vs Random", dist_linear_random),
        ("Oscillating vs Random", dist_osc_random),
    ];

    let closest = distances
        .iter()
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
        .unwrap();
    println!(
        "Most similar patterns: {} (distance: {:.6})",
        closest.0, closest.1
    );

    println!("\n=== Analysis Complete ===");
}

fn create_linear_series() -> ndarray::Array2<f64> {
    let n = 20;
    let mut data = Vec::new();

    for i in 0..n {
        let t = i as f64 * 0.1;
        data.push([t, 0.5 * t + 0.1 * t * t]); // Slightly curved line
    }

    ndarray::Array2::from_shape_vec((n, 2), data.into_iter().flatten().collect()).unwrap()
}

fn create_oscillating_series() -> ndarray::Array2<f64> {
    let n = 20;
    let mut data = Vec::new();

    for i in 0..n {
        let t = i as f64 * 0.3;
        data.push([t, (t * 2.0).sin() * 0.5]);
    }

    ndarray::Array2::from_shape_vec((n, 2), data.into_iter().flatten().collect()).unwrap()
}

fn create_random_walk() -> ndarray::Array2<f64> {
    let n = 20;
    let mut data = Vec::new();
    let mut x = 0.0;
    let mut y = 0.0;

    // Simple deterministic "random" walk
    let steps = [0.1, -0.05, 0.15, -0.1, 0.08, -0.12, 0.2, -0.15];

    for i in 0..n {
        let step_x = steps[i % steps.len()];
        let step_y = steps[(i + 3) % steps.len()];

        x += step_x;
        y += step_y;
        data.push([x, y]);
    }

    ndarray::Array2::from_shape_vec((n, 2), data.into_iter().flatten().collect()).unwrap()
}

fn analyze_signature(name: &str, sig: &LogSignature<ENotation, NotNan<f64>>) {
    println!("--- {} ---", name);

    // Compute signature norm (sum of absolute values)
    let norm: f64 = sig
        .series
        .coefficients
        .iter()
        .map(|c| c.abs().to_f64().unwrap_or(0.0))
        .sum();

    println!("  Signature norm: {:.6}", norm);

    // Find dominant terms
    let mut indexed_coeffs: Vec<(usize, f64)> = sig
        .series
        .coefficients
        .iter()
        .enumerate()
        .map(|(i, c)| (i, c.abs().to_f64().unwrap_or(0.0)))
        .collect();

    indexed_coeffs.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

    println!("  Top 3 terms:");
    for (i, (idx, val)) in indexed_coeffs.iter().take(3).enumerate() {
        if *val > 1e-10 {
            println!("    {}: {} = {:.6}", i + 1, sig.series.basis[*idx], val);
        }
    }
    println!();
}

fn compute_l2_distance(
    sig1: &LogSignature<ENotation, NotNan<f64>>,
    sig2: &LogSignature<ENotation, NotNan<f64>>,
) -> f64 {
    sig1.series
        .coefficients
        .iter()
        .zip(sig2.series.coefficients.iter())
        .map(|(a, b)| {
            let diff = **a - **b;
            diff * diff
        })
        .sum::<f64>()
        .sqrt()
}
