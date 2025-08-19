# lie-rs

A Rust library for Lie series and Baker-Campbell-Hausdorff computations.

## Overview

- Baker-Campbell-Hausdorff (BCH) series computation
- Lie series representation and manipulation
- Generic numeric types (rationals, floats, integers)

## Quick Start

Add this to your `Cargo.toml`:

```toml
[dependencies]
lie-rs = "0.1.0"
```

### Basic BCH Computation

```rust
use lie_rs::prelude::*;
use lyndon_rs::{generators::ENotation, prelude::*};

// Create a BCH series generator for 2 generators up to degree 5
let basis = LyndonBasis::<ENotation>::new(2, Sort::Lexicographical);
let generator = BchSeriesGenerator::<ENotation>::new(basis, 5);

// Generate the BCH series
let bch_series = generator.generate_lie_series();

println!("BCH series with {} terms", bch_series.len());

// Print the first few terms
for (word, coeff) in bch_series.iter().take(10) {
    println!("  {:?}: {}", word, coeff);
}
```

### Working with Rational Arithmetic

```rust
use lie_rs::prelude::*;
use lyndon_rs::{generators::ENotation, prelude::*};
use num_rational::Ratio;

// Use exact rational arithmetic for precise computation
let basis = LyndonBasis::<ENotation>::new(3, Sort::Lexicographical);
let generator = BchSeriesGenerator::<ENotation>::new(basis, 4);
let bch_series: LieSeries<ENotation, Ratio<i64>> = generator.generate_lie_series();

// All coefficients are exact rational numbers
for (word, coeff) in bch_series.iter() {
    println!("{:?}: {}/{}", word, coeff.numer(), coeff.denom());
}
```

### Creating and Manipulating Lie Series

```rust
use lie_rs::prelude::*;
use lyndon_rs::prelude::*;

// Create a custom Lie series with basis and coefficients
let word1 = LyndonWord::try_from(vec![ENotation::new(1)]).unwrap();
let word2 = LyndonWord::try_from(vec![ENotation::new(2)]).unwrap();
let word12 = LyndonWord::try_from(vec![ENotation::new(1), ENotation::new(2)]).unwrap();

let basis = vec![word1, word2, word12];
let coefficients = vec![1.0, 1.0, 0.5];

let series = LieSeries::<ENotation, f64>::new(basis, coefficients);

println!("Custom series: {:?}", series);

// Series support algebraic operations
let doubled = &series + &series;  // Addition
println!("Doubled series: {:?}", doubled);
```

### Computing Commutators in Lie Series

```rust
use lie_rs::prelude::*;
use lyndon_rs::{generators::ENotation, prelude::*};
use commutator_rs::Commutator;

// Create two simple Lie series
let x = LyndonWord::try_from(vec![ENotation::new(1)]).unwrap();
let y = LyndonWord::try_from(vec![ENotation::new(2)]).unwrap();

let x_series = LieSeries::<ENotation, f64>::new(vec![x.clone()], vec![1.0]);
let y_series = LieSeries::<ENotation, f64>::new(vec![y.clone()], vec![1.0]);

// Compute the commutator [x, y]
let commutator_series = x_series.commutator(&y_series);
println!("Commutator [x, y]: {:?}", commutator_series);
```

## Key Types

- `LieSeriesGenerator`: Trait for generating Lie series from various inputs
- `BchSeriesGenerator`: Specialized generator for Baker-Campbell-Hausdorff computations
- `LieSeries<T, U>`: Generic container for Lie algebraic series
- `RootedTree<T>`: Trees used for BCH computations
