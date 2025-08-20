//! Formal indeterminate operations example.
//!
//! This example demonstrates working with formal indeterminates,
//! which represent products of symbols with scalar coefficients.

use commutator_rs::prelude::*;
use std::collections::HashMap;

fn main() {
    println!("=== Formal Indeterminate Operations ===\n");

    // Basic formal indeterminate creation
    demonstrate_basic_indeterminates();

    // Operations on formal indeterminates
    println!("\n=== Indeterminate Operations ===");
    demonstrate_indeterminate_operations();

    // Converting commutators to formal indeterminates
    println!("\n=== Commutator to Indeterminate Conversion ===");
    demonstrate_commutator_conversion();

    // Building algebraic expressions
    println!("\n=== Algebraic Expressions ===");
    demonstrate_algebraic_expressions();

    println!("\n=== Example Complete ===");
}

fn demonstrate_basic_indeterminates() {
    println!("Creating and displaying formal indeterminates:");
    println!();

    // Single symbol indeterminates
    let x = FormalIndeterminate::new(vec!["x"], 1.0);
    let y = FormalIndeterminate::new(vec!["y"], 2.0);

    println!("Single symbols:");
    println!("  x with coefficient 1.0: {}", x);
    println!("  y with coefficient 2.0: {}", y);
    println!();

    // Multi-symbol indeterminates (products)
    let xy = FormalIndeterminate::new(vec!["x", "y"], 1.0);
    let xyz = FormalIndeterminate::new(vec!["x", "y", "z"], -0.5);

    println!("Product symbols:");
    println!("  xy with coefficient 1.0: {}", xy);
    println!("  xyz with coefficient -0.5: {}", xyz);
    println!();

    // With different coefficient types
    let rational_term = FormalIndeterminate::new(vec!['a', 'b'], 3);
    println!("With integer coefficients:");
    println!("  ab with coefficient 3: {}", rational_term);
}

fn demonstrate_indeterminate_operations() {
    println!("Operations on formal indeterminates:");
    println!();

    let x = FormalIndeterminate::new(vec!["x"], 2.0);
    let y = FormalIndeterminate::new(vec!["y"], 3.0);

    println!("Starting indeterminates:");
    println!("  x: {}", x);
    println!("  y: {}", y);
    println!();

    // Multiplication combines symbols and multiplies coefficients
    let product = &x * &y;
    println!("Multiplication x * y:");
    println!("  Result: {}", product);
    println!("  Combined symbols: {:?}", product.symbols);
    println!("  Combined coefficient: {}", product.coefficient);
    println!();

    // Scalar multiplication
    let scaled = &x * 5.0;
    println!("Scalar multiplication x * 5.0:");
    println!("  Result: {}", scaled);
    println!();

    // Negation
    let negated = -&x;
    println!("Negation -x:");
    println!("  Result: {}", negated);
    println!();

    // More complex products
    let z = FormalIndeterminate::new(vec!["z"], 1.0);
    let triple_product = &(&x * &y) * &z;
    println!("Triple product (x * y) * z:");
    println!("  Result: {}", triple_product);
    println!("  Symbols: {:?}", triple_product.symbols);
}

fn demonstrate_commutator_conversion() {
    println!("Converting commutator expressions to formal indeterminates:");
    println!();

    // Create a simple commutator
    let a = CommutatorTerm::Atom {
        coefficient: 1,
        atom: 'A',
    };
    let b = CommutatorTerm::Atom {
        coefficient: 1,
        atom: 'B',
    };

    let comm_ab = a.commutator(&b);
    println!("Commutator [A, B]:");
    println!("  Structure: {:?}", comm_ab);
    println!();

    // Convert to formal indeterminates
    let indeterminates = Vec::<FormalIndeterminate<char, i32>>::from(&comm_ab);

    println!("Converted to formal indeterminates:");
    for (i, indeterminate) in indeterminates.iter().enumerate() {
        println!("  Term {}: {}", i + 1, indeterminate);
        println!("    Coefficient: {}", indeterminate.coefficient);
        println!("    Symbols: {:?}", indeterminate.symbols);
    }
    println!();

    // Show the expansion: [A,B] = AB - BA
    println!("This represents the expansion [A,B] = AB - BA");

    // Try a more complex commutator
    println!("\nMore complex example:");
    let c = CommutatorTerm::Atom {
        coefficient: 1,
        atom: 'C',
    };
    let nested_comm = comm_ab.commutator(&c);

    let nested_indeterminates = Vec::<FormalIndeterminate<char, i32>>::from(&nested_comm);
    println!(
        "[[A,B], C] expands to {} terms:",
        nested_indeterminates.len()
    );

    for (i, indeterminate) in nested_indeterminates.iter().enumerate().take(8) {
        println!("  Term {}: {}", i + 1, indeterminate);
    }
    if nested_indeterminates.len() > 8 {
        println!("  ... and {} more terms", nested_indeterminates.len() - 8);
    }
}

fn demonstrate_algebraic_expressions() {
    println!("Building complex algebraic expressions:");
    println!();

    // Create a collection of terms representing a polynomial-like expression
    let mut expression = HashMap::new();

    // Add various monomial terms
    let x = FormalIndeterminate::new(vec!["x"], 1.0);
    let y = FormalIndeterminate::new(vec!["y"], 1.0);
    let x2 = &x * &x; // x²
    let xy = &x * &y; // xy
    let y2 = &y * &y; // y²

    // Store with string keys for display
    expression.insert("x".to_string(), x);
    expression.insert("y".to_string(), y);
    expression.insert("x²".to_string(), x2);
    expression.insert("xy".to_string(), xy);
    expression.insert("y²".to_string(), y2);

    println!("Algebraic expression terms:");
    for (name, term) in &expression {
        println!("  {}: {} (symbols: {:?})", name, term, term.symbols);
    }
    println!();

    // Demonstrate with different coefficient types
    println!("Working with fractional coefficients:");

    let fractional_x = FormalIndeterminate::new(vec!["x"], 0.5);
    let fractional_y = FormalIndeterminate::new(vec!["y"], 0.75);

    println!("  x with coefficient 0.5: {}", fractional_x);
    println!("  y with coefficient 0.75: {}", fractional_y);

    let fractional_product = &fractional_x * &fractional_y;
    println!("  Product xy: {}", fractional_product);
    println!("  Product coefficient: {}", fractional_product.coefficient);

    // Show arithmetic computation
    println!();
    println!("Arithmetic computation:");
    println!("  0.5 * 0.75 = {}", fractional_product.coefficient);

    // Complex symbolic computation example
    println!("\nSymbolic computation example:");
    demonstrate_symbolic_computation();
}

fn demonstrate_symbolic_computation() {
    // Create generators for a small algebra
    let e1 = FormalIndeterminate::new(vec![1], 1);
    let e2 = FormalIndeterminate::new(vec![2], 1);
    let e3 = FormalIndeterminate::new(vec![3], 1);

    println!("Generators:");
    println!("  e₁: {}", e1);
    println!("  e₂: {}", e2);
    println!("  e₃: {}", e3);
    println!();

    // Build various products
    let products = vec![
        ("e₁e₂", &e1 * &e2),
        ("e₁e₃", &e1 * &e3),
        ("e₂e₃", &e2 * &e3),
        ("e₁e₂e₃", &(&e1 * &e2) * &e3),
        ("e₂e₁e₃", &(&e2 * &e1) * &e3),
    ];

    println!("Various products:");
    for (name, product) in products {
        println!(
            "  {}: {} (length: {})",
            name,
            product,
            product.symbols.len()
        );
    }
    println!();

    // Show how this relates to non-commutative algebra
    println!("In non-commutative algebra:");
    println!("  e₁e₂ ≠ e₂e₁ (generally)");
    println!("  Order of symbols matters");
    println!("  Formal indeterminates track the exact order");

    // Demonstrate coefficient accumulation
    println!("\nCoefficient calculations:");
    let term1 = FormalIndeterminate::new(vec![1, 2], 2);
    let term2 = FormalIndeterminate::new(vec![1, 2], 3);
    // In a real algebra, we'd add these: (2 + 3) * e₁e₂ = 5 * e₁e₂

    println!("  2·e₁e₂: {}", term1);
    println!("  3·e₁e₂: {}", term2);
    println!("  Sum would be: 5·e₁e₂ (requires addition implementation)");
}
