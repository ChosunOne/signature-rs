# commutator-rs

A Rust library for commutator operations.

## Overview

- Commutator operations: `[A, B] = AB - BA`
- Commutator terms for algebraic expressions
- Formal indeterminates for symbolic computation
- Expression trees for nested commutators

## Quick Start

Add this to your `Cargo.toml`:

```toml
[dependencies]
commutator-rs = "0.1.0"
```

### Basic Commutator Operations

```rust
use commutator_rs::prelude::*;

// Create formal indeterminates
let x = FormalIndeterminate::new(vec!["x"], 1.0);
let y = FormalIndeterminate::new(vec!["y"], 1.0);

// Compute commutator [x, y] using the trait
let result = x.commutator(&y);
println!("Commutator [x, y]: {:?}", result);

// For simple types, commutator of identical elements is zero
let zero_result = 5.commutator(&5);
println!("Commutator [5, 5]: {:?}", zero_result); // Should be 0
```

### Using the `comm!` Macro

For convenient syntax, use the `comm!` macro instead of calling `.commutator()` directly:

```rust
use commutator_rs::{CommutatorTerm, Commutator, comm};

let x = CommutatorTerm::Atom { coefficient: 1, atom: 'x' };
let y = CommutatorTerm::Atom { coefficient: 1, atom: 'y' };

// These are equivalent:
let result1 = x.commutator(&y);
let result2 = comm![x, y];

// Nested commutators are easy to read:
let a = CommutatorTerm::<i32, char>::from('a');
let b = CommutatorTerm::<i32, char>::from('b'); 
let c = CommutatorTerm::<i32, char>::from('c');

let nested = comm![comm![a, b], c]; // [[a, b], c]
```

### Working with Commutator Terms

```rust
use commutator_rs::prelude::*;
use lyndon_rs::generators::ENotation;

// Create commutator terms from generators
let e1 = ENotation::new(1);
let e2 = ENotation::new(2);
let e3 = ENotation::new(3);

// Build a simple commutator term
let term1 = CommutatorTerm::Atom { coefficient: 1, atom: e1.clone() };
let term2 = CommutatorTerm::Atom { coefficient: 1, atom: e2.clone() };

// Create nested commutator: [[e1, e2], e3]
let inner_comm = CommutatorTerm::Expression {
    coefficient: 1,
    left: Box::new(term1),
    right: Box::new(term2),
};

let term3 = CommutatorTerm::Atom { coefficient: 1, atom: e3 };
let outer_comm = CommutatorTerm::Expression {
    coefficient: 1,
    left: Box::new(inner_comm),
    right: Box::new(term3),
};

println!("Nested commutator: {:?}", outer_comm);
```

### Creating Complex Expressions

```rust
use commutator_rs::prelude::*;
use std::collections::HashMap;

// Create a more complex example with coefficients
let mut expression = HashMap::new();

// Create formal indeterminates
let x = FormalIndeterminate::new("x", 1.0);
let y = FormalIndeterminate::new("y", 2.0);
let z = FormalIndeterminate::new("z", 0.5);

// Compute various commutators
let xy_comm = x.commutator(&y);
let xz_comm = x.commutator(&z);
let yz_comm = y.commutator(&z);

println!("[x, y] = {:?}", xy_comm);
println!("[x, z] = {:?}", xz_comm);
println!("[y, z] = {:?}", yz_comm);

// These can be used in further algebraic computations
let nested = xy_comm.commutator(&z);
println!("[[x, y], z] = {:?}", nested);
```

### Integration with Lyndon Words

```rust
use commutator_rs::prelude::*;
use lyndon_rs::prelude::*;

// Create commutator terms using generators from Lyndon words
let generators = vec![ENotation::new(1), ENotation::new(2)];

// Try to create a Lyndon word 
if let Ok(lyndon_word) = LyndonWord::try_from(generators) {
    println!("Valid Lyndon word: {:?}", lyndon_word);
    
    // Create commutator terms manually from the Lyndon word's letters
    let term1 = CommutatorTerm::Atom { coefficient: 1, atom: lyndon_word.letters[0].clone() };
    let term2 = CommutatorTerm::Atom { coefficient: 1, atom: lyndon_word.letters[1].clone() };
    let comm_term = term1.commutator(&term2);
    println!("Commutator from Lyndon word: {:?}", comm_term);
}
```

### Custom Types with Commutators

```rust
use commutator_rs::Commutator;
use std::ops::{Mul, Sub};

// Define a custom matrix-like type
#[derive(Clone, Debug, PartialEq)]
struct Matrix2x2 {
    a: f64, b: f64,
    c: f64, d: f64,
}

impl Matrix2x2 {
    fn new(a: f64, b: f64, c: f64, d: f64) -> Self {
        Self { a, b, c, d }
    }
}

impl Mul for Matrix2x2 {
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

impl Sub for Matrix2x2 {
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

impl Commutator for Matrix2x2 {}

// Usage
let m1 = Matrix2x2::new(1.0, 2.0, 3.0, 4.0);
let m2 = Matrix2x2::new(0.0, 1.0, -1.0, 0.0);

let comm = m1.commutator(&m2);
println!("Matrix commutator: {:?}", comm);
```

## Key Types

- `Commutator`: Trait defining the commutator operation
- `CommutatorTerm<T, U>`: Represents either atomic terms or nested commutator expressions
- `FormalIndeterminate<T, U>`: Abstract algebraic elements for symbolic computation
- `comm!`: Macro for commutator syntax
