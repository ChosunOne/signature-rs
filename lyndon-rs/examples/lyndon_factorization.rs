//! Lyndon word factorization and analysis example.
//!
//! This example demonstrates advanced features of Lyndon words including
//! factorization and their role in free Lie algebra bases.

use lyndon_rs::prelude::*;
use std::collections::HashMap;

fn main() {
    println!("=== Lyndon Word Factorization Example ===\n");

    // Create basis for demonstration
    let basis = LyndonBasis::<ENotation>::new(2, Sort::Lexicographical);
    let words = basis.generate_basis(6);

    println!("Generated {} Lyndon words up to length 6", words.len());
    println!("Using binary alphabet (size 2)\n");

    // Group words by length
    let mut words_by_length: HashMap<usize, Vec<&LyndonWord<ENotation>>> = HashMap::new();
    for word in &words {
        words_by_length
            .entry(word.letters.len())
            .or_default()
            .push(word);
    }

    println!("=== Words by Length ===");
    for length in 1..=6 {
        if let Some(words_of_length) = words_by_length.get(&length) {
            println!("Length {}: {} words", length, words_of_length.len());
            for (i, word) in words_of_length.iter().enumerate() {
                println!("  {}: {}", i + 1, word);
            }
            println!();
        }
    }

    println!("=== Lyndon Word Construction Analysis ===");
    analyze_word_construction(&words);

    println!("\n=== Free Lie Algebra Connection ===");
    demonstrate_lie_algebra_connection(&basis);

    println!("\n=== Möbius Function Application ===");
    demonstrate_moebius_function();

    println!("\n=== Example Complete ===");
}

fn analyze_word_construction(words: &[LyndonWord<ENotation>]) {
    println!("Analyzing how Lyndon words are constructed:");
    println!();

    // Look at words and their structure
    for word in words.iter().take(12) {
        if word.letters.len() > 1 {
            println!("Word: {}", word);
            println!("  Letters: {:?}", word.letters);

            // Try to find factorization as (u,v) where word = uv and u,v are Lyndon
            analyze_factorizations(word, words);
            println!();
        }
    }
}

fn analyze_factorizations(word: &LyndonWord<ENotation>, all_words: &[LyndonWord<ENotation>]) {
    let n = word.letters.len();
    let mut factorizations = Vec::new();

    // Try all possible splits
    for i in 1..n {
        let left_part = &word.letters[0..i];
        let right_part = &word.letters[i..n];

        // Check if both parts are Lyndon words
        let left_is_lyndon = all_words.iter().any(|w| w.letters == left_part);
        let right_is_lyndon = all_words.iter().any(|w| w.letters == right_part);

        if left_is_lyndon && right_is_lyndon {
            factorizations.push((left_part.to_vec(), right_part.to_vec()));
        }
    }

    if factorizations.is_empty() {
        println!("  No Lyndon factorizations found (likely primitive)");
    } else {
        println!("  Lyndon factorizations:");
        for (i, (left, right)) in factorizations.iter().enumerate() {
            println!("    {}: {:?} + {:?}", i + 1, left, right);
        }
    }
}

fn demonstrate_lie_algebra_connection(basis: &LyndonBasis<ENotation>) {
    println!("Connection to Free Lie Algebra:");
    println!();

    // The number of Lyndon words of each degree corresponds to
    // the dimensions of graded components of the free Lie algebra

    let max_degree = 8;
    let counts = basis.number_of_words_per_degree(max_degree);

    println!("Lyndon word counts = Free Lie algebra dimensions:");
    println!(
        "(For free Lie algebra on {} generators)",
        basis.alphabet_size
    );
    println!();

    let mut total = 0;
    for (degree, count) in counts.iter().enumerate().skip(1) {
        if *count > 0 {
            total += count;
            println!(
                "  Degree {}: {} words → dim(L_{}) = {}",
                degree, count, degree, count
            );
        }
    }

    println!("\nTotal dimension up to degree {}: {}", max_degree, total);

    // Show the connection to the Möbius function
    println!("\nThis count is given by the formula:");
    println!("  |Lyndon_n(A)| = (1/n) * Σ μ(d) * |A|^(n/d)");
    println!("where the sum is over divisors d of n, and μ is the Möbius function");
}

fn demonstrate_moebius_function() {
    println!("Möbius function values for small integers:");
    println!();

    for n in 1..=12 {
        let mu = lyndon_rs::lyndon::moebius_mu(n);
        let interpretation = match mu[0] {
            // Get the first element
            1 => "square-free with even number of prime factors",
            -1 => "square-free with odd number of prime factors",
            0 => "not square-free (has squared prime factor)",
            _ => "unexpected value",
        };

        println!("  μ({:2}) = {:2}  ({})", n, mu[0], interpretation);
    }

    println!();
    println!("The Möbius function is used in the formula for counting Lyndon words:");

    // Demonstrate the counting formula for a few cases
    let alphabet_size = 2;
    println!("\nFor alphabet size {}:", alphabet_size);

    for n in 1..=6 {
        let divisors = get_divisors(n);
        let mut sum = 0;
        let mut formula_parts = Vec::new();

        for d in divisors {
            let mu_d = lyndon_rs::lyndon::moebius_mu(d);
            let term = mu_d[0] * (alphabet_size as i64).pow((n / d) as u32);
            sum += term;

            if term != 0 {
                let sign = if term > 0 { "+" } else { "" };
                formula_parts.push(format!("{}μ({})×{}^{}", sign, d, alphabet_size, n / d));
            }
        }

        let count = sum / n as i64;
        let formula = formula_parts.join(" ");

        println!("  n={}: ({}) / {} = {}", n, formula, n, count);
    }
}

fn get_divisors(n: usize) -> Vec<usize> {
    let mut divisors = Vec::new();
    for i in 1..=n {
        if n % i == 0 {
            divisors.push(i);
        }
    }
    divisors
}
