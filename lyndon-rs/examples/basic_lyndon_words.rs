//! Basic Lyndon word generation example.
//!
//! This example demonstrates how to generate Lyndon words using the lyndon-rs library
//! and explore their properties.

use lyndon_rs::prelude::*;

fn main() {
    println!("=== Lyndon Words Generation Example ===\n");

    // Create a Lyndon basis with different configurations
    let basis_small = LyndonBasis::<ENotation>::new(2, Sort::Lexicographical);
    let basis_large = LyndonBasis::<ENotation>::new(3, Sort::Lexicographical);

    println!("=== Small Alphabet (size 2) ===");
    demonstrate_basis(&basis_small, "Binary", 5);

    println!("\n=== Larger Alphabet (size 3) ===");
    demonstrate_basis(&basis_large, "Ternary", 4);

    println!("\n=== Sorting Methods Comparison ===");
    compare_sorting_methods();

    println!("\n=== Lyndon Word Properties ===");
    explore_lyndon_properties();

    println!("\n=== Example Complete ===");
}

fn demonstrate_basis(basis: &LyndonBasis<ENotation>, name: &str, max_length: usize) {
    println!("--- {} Alphabet (size {}) ---", name, basis.alphabet_size);

    // Count words by degree
    let counts = basis.number_of_words_per_degree(max_length);
    println!("Number of Lyndon words by degree:");
    for (degree, count) in counts.iter().enumerate().skip(1) {
        if *count > 0 {
            println!("  Degree {}: {} words", degree, count);
        }
    }

    let total: usize = counts.iter().sum();
    println!("Total words up to degree {}: {}", max_length, total);
    println!();

    // Generate and display words
    let words = basis.generate_basis(max_length);

    println!("Generated Lyndon words:");
    for (i, word) in words.iter().enumerate() {
        if i < 15 {
            // Show first 15 words
            println!("  {}: {}", i + 1, word);
        }
    }

    if words.len() > 15 {
        println!("  ... and {} more words", words.len() - 15);
    }
    println!();
}

fn compare_sorting_methods() {
    let alphabet_size = 2;
    let max_length = 4;

    let basis_lex = LyndonBasis::<ENotation>::new(alphabet_size, Sort::Lexicographical);
    let basis_topo = LyndonBasis::<ENotation>::new(alphabet_size, Sort::Topological);

    let words_lex = basis_lex.generate_basis(max_length);
    let words_topo = basis_topo.generate_basis(max_length);

    println!(
        "Comparing sorting methods (alphabet size {}, max length {}):",
        alphabet_size, max_length
    );
    println!();

    println!("Lexicographical ordering:");
    for (i, word) in words_lex.iter().take(10).enumerate() {
        println!("  {}: {}", i + 1, word);
    }

    println!("\nTopological ordering:");
    for (i, word) in words_topo.iter().take(10).enumerate() {
        println!("  {}: {}", i + 1, word);
    }

    // Verify they contain the same words (just different order)
    let same_content =
        words_lex.len() == words_topo.len() && words_lex.iter().all(|w| words_topo.contains(w));

    println!(
        "\nSame word content: {}",
        if same_content { "✓ Yes" } else { "✗ No" }
    );
}

fn explore_lyndon_properties() {
    println!("=== Lyndon Word Properties ===");

    // Create some example words
    let basis = LyndonBasis::<ENotation>::new(3, Sort::Lexicographical);
    let words = basis.generate_basis(3);

    println!("Exploring properties of first few words:");

    for (i, word) in words.iter().take(8).enumerate() {
        println!("\nWord {}: {}", i + 1, word);
        println!("  Length: {}", word.letters.len());
        println!("  Letters: {:?}", word.letters);

        // Check if it's primitive (not a power of a shorter word)
        let is_primitive = check_if_primitive(word);
        println!(
            "  Primitive: {}",
            if is_primitive { "✓ Yes" } else { "✗ No" }
        );

        // Show lexicographic comparison with rotations
        demonstrate_lexicographic_property(word);
    }
}

fn check_if_primitive(word: &LyndonWord<ENotation>) -> bool {
    let n = word.letters.len();
    if n <= 1 {
        return true;
    }

    // Check if word is a repetition of a shorter period
    for period in 1..n {
        if n % period == 0 {
            let _repetitions = n / period;
            let mut is_repetition = true;

            for i in 0..n {
                if word.letters[i] != word.letters[i % period] {
                    is_repetition = false;
                    break;
                }
            }

            if is_repetition {
                return false;
            }
        }
    }

    true
}

fn demonstrate_lexicographic_property(word: &LyndonWord<ENotation>) {
    if word.letters.len() <= 1 {
        println!("    (Too short for rotation comparison)");
        return;
    }

    // Generate all rotations
    let n = word.letters.len();
    let mut rotations = Vec::new();

    for i in 0..n {
        let mut rotation = Vec::new();
        for j in 0..n {
            rotation.push(word.letters[(i + j) % n].clone());
        }
        rotations.push(rotation);
    }

    // Check if original is lexicographically smallest
    let original = &word.letters;
    let is_minimal = rotations.iter().all(|rot| original <= rot);

    println!(
        "    Minimal rotation: {}",
        if is_minimal { "✓ Yes" } else { "✗ No" }
    );

    // Show first few rotations for illustration
    if rotations.len() <= 4 {
        println!("    All rotations:");
        for (i, rot) in rotations.iter().enumerate() {
            let marker = if rot == original { " ← original" } else { "" };
            println!("      {}: {:?}{}", i + 1, rot, marker);
        }
    }
}
