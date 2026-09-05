//! Basic commutator operations example.
//!
//! This example demonstrates fundamental commutator operations using
//! the commutator-rs library, including the convenient comm! macro.

use commutator_rs::comm;
use commutator_rs::prelude::*;

/// A minimal 2x2 matrix implementing `Mul`/`Add`/`Sub`, so the blanket
/// `Commutator` impl applies and Jacobi-identity sums can be computed
/// numerically.
#[derive(Clone, Debug, PartialEq)]
struct Matrix2x2 {
    a: f64,
    b: f64,
    c: f64,
    d: f64,
}

impl Matrix2x2 {
    fn new(a: f64, b: f64, c: f64, d: f64) -> Self {
        Self { a, b, c, d }
    }

    fn show(&self, label: &str) {
        println!("  {label}: [{:.1} {:.1}]", self.a, self.b);
        println!("        [{:.1} {:.1}]", self.c, self.d);
    }
}

impl std::ops::Mul for Matrix2x2 {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        Self {
            a: self.a * other.a + self.b * other.c,
            b: self.a * other.b + self.b * other.d,
            c: self.c * other.a + self.d * other.c,
            d: self.c * other.b + self.d * other.d,
        }
    }
}

impl std::ops::Add for Matrix2x2 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            a: self.a + other.a,
            b: self.b + other.b,
            c: self.c + other.c,
            d: self.d + other.d,
        }
    }
}

impl std::ops::Sub for Matrix2x2 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            a: self.a - other.a,
            b: self.b - other.b,
            c: self.c - other.c,
            d: self.d - other.d,
        }
    }
}

fn main() {
    println!("=== Basic Commutator Operations ===\n");

    // Demonstrate with CommutatorTerm
    demonstrate_commutator_terms();

    // Demonstrate with numeric types
    println!("\n=== Numeric Commutators ===");
    demonstrate_numeric_commutators();

    // Demonstrate the comm! macro
    println!("\n=== Using the comm! Macro ===");
    demonstrate_comm_macro();

    // Show commutator properties
    println!("\n=== Commutator Properties ===");
    demonstrate_commutator_properties();

    // Complex nested commutators
    println!("\n=== Nested Commutators ===");
    demonstrate_nested_commutators();

    println!("\n=== Example Complete ===");
}

fn demonstrate_commutator_terms() {
    println!("Creating and computing CommutatorTerm commutators:");
    println!();

    // Create basic atomic terms
    let x = CommutatorTerm::Atom {
        coefficient: 1,
        atom: 'x',
    };
    let y = CommutatorTerm::Atom {
        coefficient: 1,
        atom: 'y',
    };
    let z = CommutatorTerm::Atom {
        coefficient: 1,
        atom: 'z',
    };

    println!("Created atoms:");
    println!("  x = {x}");
    println!("  y = {y}");
    println!("  z = {z}");
    println!();

    // Compute basic commutators
    let xy_comm = x.commutator(&y);
    let xz_comm = x.commutator(&z);
    let yz_comm = y.commutator(&z);

    println!("Basic commutators:");
    println!("  [x, y] = {xy_comm}");
    println!("  [x, z] = {xz_comm}");
    println!("  [y, z] = {yz_comm}");
    println!();

    // Show commutator structure
    if let CommutatorTerm::Expression {
        coefficient,
        left,
        right,
        degree,
    } = &xy_comm
    {
        println!("Structure of [x, y]:");
        println!("  Coefficient: {coefficient}");
        println!("  Left: {left}");
        println!("  Right: {right}");
        println!("  Degree: {degree}");
    }
}

fn demonstrate_numeric_commutators() {
    // For numeric types, [a,b] = ab - ba
    println!("Numeric commutators (matrices, etc.):");
    println!();

    // Demonstrate with simple numbers (commutative, so result is 0)
    let a = 3;
    let b = 5;
    let comm_result = a.commutator(&b);
    println!("  [3, 5] = 3×5 - 5×3 = {comm_result} (numbers commute)");

    // Demonstrate concept with a custom matrix-like structure
    demonstrate_matrix_commutators();
}

fn demonstrate_matrix_commutators() {
    println!("\nMatrix-like commutators ([A,B] = AB - BA):\n");

    // Matrix2x2 automatically gets the Commutator trait from the blanket
    // implementation over Mul + Sub types.
    let m1 = Matrix2x2::new(1.0, 2.0, 3.0, 4.0);
    let m2 = Matrix2x2::new(0.0, 1.0, -1.0, 0.0);

    m1.show("Matrix 1");
    m2.show("Matrix 2");

    let comm = m1.commutator(&m2);
    comm.show("[M1, M2] ");
}

fn demonstrate_comm_macro() {
    println!("The comm! macro provides convenient syntax:");
    println!();

    let x = CommutatorTerm::<i32, char>::from('x');
    let y = CommutatorTerm::<i32, char>::from('y');
    let z = CommutatorTerm::<i32, char>::from('z');

    // Compare syntax
    let method_result = x.commutator(&y);
    let macro_result = comm![x, y];

    println!("  Method syntax: x.commutator(&y)");
    println!("  Macro syntax:  comm![x, y]");
    println!("  Results equal: {}", method_result == macro_result);
    println!();

    // Nested commutators are much cleaner with macro
    let nested_method = (x.commutator(&y)).commutator(&z);
    let nested_macro = comm![comm![x, y], z];

    println!("  Nested method: (x.commutator(&y)).commutator(&z)");
    println!("  Nested macro:  comm![comm![x, y], z]");
    println!("  Results equal: {}", nested_method == nested_macro);
}

fn demonstrate_commutator_properties() {
    println!("Fundamental commutator properties:");
    println!();

    let x = CommutatorTerm::<i32, char>::from('x');
    let y = CommutatorTerm::<i32, char>::from('y');
    let z = CommutatorTerm::<i32, char>::from('z');

    // 1. Anti-commutativity: [x,y] = -[y,x]
    let xy = comm![x, y];
    let yx = comm![y, x];
    let neg_yx = -yx;

    println!("1. Anti-commutativity: [x,y] = -[y,x]");
    println!("   [x,y]  = {xy}");
    println!("   [y,x]  = {}", comm![y, x]);
    println!("  -[y,x]  = {neg_yx}");
    println!("   Equal: {}", xy == neg_yx);
    println!();

    // 2. Self-commutator is zero: [x,x] = 0
    let xx = comm![x, x];
    println!("2. Self-commutator: [x,x] = 0");
    println!("   [x,x] = {xx}");
    println!("   Is zero: {}", xx.is_zero());
    println!();

    // 3. Jacobi identity: [x,[y,z]] + [y,[z,x]] + [z,[x,y]] = 0.
    // CommutatorTerm has no Add, so the identity is verified numerically
    // with matrices (where sums exist) via the same commutator trait.
    println!("3. Jacobi identity: [x,[y,z]] + [y,[z,x]] + [z,[x,y]] = 0");
    let term1 = comm![x, comm![y, z]];
    let term2 = comm![y, comm![z, x]];
    let term3 = comm![z, comm![x, y]];

    println!("   [x,[y,z]] = {term1}");
    println!("   [y,[z,x]] = {term2}");
    println!("   [z,[x,y]] = {term3}");

    let mx = Matrix2x2::new(0.0, 1.0, -1.0, 0.0);
    let my = Matrix2x2::new(0.0, 0.0, 1.0, 0.0);
    let mz = Matrix2x2::new(1.0, 0.0, 0.0, -1.0);
    let jacobi = mx
        .commutator(&my.commutator(&mz))
        + my.commutator(&mz.commutator(&mx))
        + mz.commutator(&mx.commutator(&my));
    let is_zero = jacobi.a.abs() < 1e-12
        && jacobi.b.abs() < 1e-12
        && jacobi.c.abs() < 1e-12
        && jacobi.d.abs() < 1e-12;
    println!("   Numeric check (2x2 matrices): sum is zero: {is_zero}");
}

fn demonstrate_nested_commutators() {
    println!("Working with complex nested commutator expressions:");
    println!();

    let a = CommutatorTerm::<i32, char>::from('a');
    let b = CommutatorTerm::<i32, char>::from('b');
    let c = CommutatorTerm::<i32, char>::from('c');

    // Build increasingly complex expressions
    println!("Building nested commutators:");

    let ab = comm![a, b];
    println!("  [a,b] has degree: {}", ab.degree());

    let ab_c = comm![ab, c];
    println!("  [[a,b],c] has degree: {}", ab_c.degree());

    let bc = comm![b, c];
    let a_bc = comm![a, bc];
    println!("  [a,[b,c]] has degree: {}", a_bc.degree());

    // Very nested example
    let complex = comm![comm![a, comm![b, c]], comm![b, a]];
    println!("  [[a,[b,c]],[b,a]] has degree: {}", complex.degree());

    println!();
    println!("The degree counts the total number of atomic elements in the expression.");

    // Show structure of a nested commutator
    println!("\nStructure of [[a,b],c]:");
    print_commutator_structure(&ab_c, 0);
}

fn print_commutator_structure(term: &CommutatorTerm<i32, char>, indent: usize) {
    let spaces = "  ".repeat(indent);
    match term {
        CommutatorTerm::Atom { coefficient, atom } => {
            println!("{}Atom: {} * '{}'", spaces, coefficient, atom);
        }
        CommutatorTerm::Expression {
            coefficient,
            left,
            right,
            ..
        } => {
            println!("{}Expression (coeff: {}):", spaces, coefficient);
            println!("{}  Left:", spaces);
            print_commutator_structure(left, indent + 2);
            println!("{}  Right:", spaces);
            print_commutator_structure(right, indent + 2);
        }
    }
}
