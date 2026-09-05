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
use ordered_float::NotNan;

// Create a BCH series generator for 2 generators up to degree 5
let basis = LyndonBasis::<ENotation>::new(2, Sort::Lexicographical);
let generator = BchSeriesGenerator::<ENotation>::new(basis, 5);

// Generate the BCH series. Coefficient types need `Eq + Hash + Ord`
// (they key the structure-constant tables): wrap raw floats in
// `NotNan`, or use integers/rationals.
let bch_series: LieSeries<ENotation, NotNan<f64>> = generator.generate_lie_series();

println!("BCH series with {} terms", bch_series.basis.len());

// Print the first few terms
for (word, coeff) in bch_series.basis.iter().zip(bch_series.coefficients.iter()).take(10) {
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
for (word, coeff) in bch_series.basis.iter().zip(bch_series.coefficients.iter()) {
    if coeff.numer() != &0 {
        println!("{:?}: {}/{}", word, coeff.numer(), coeff.denom());
    }
}
```

### Creating and Manipulating Lie Series

```rust
use lie_rs::prelude::*;
use lyndon_rs::prelude::*;
use ordered_float::NotNan;

let e = ENotation::alphabet(2);

// Create a custom Lie series with basis and coefficients
let word1 = LyndonWord::try_from(vec![e[0].clone()]).unwrap();
let word2 = LyndonWord::try_from(vec![e[1].clone()]).unwrap();
let word12 = LyndonWord::try_from(vec![e[0].clone(), e[1].clone()]).unwrap();

let basis = vec![word1, word2, word12];
// Coefficient types need `Eq + Hash + Ord`; wrap raw floats in `NotNan`.
let coefficients = vec![
    NotNan::new(1.0).unwrap(),
    NotNan::new(1.0).unwrap(),
    NotNan::new(0.5).unwrap(),
];

let series = LieSeries::<ENotation, NotNan<f64>>::new(basis, coefficients);

println!("Custom series: {:?}", series);

// Series support algebraic operations
let doubled = &series + &series;  // Addition
println!("Doubled series: {:?}", doubled);
```

### Computing Commutators in Lie Series

```rust
use commutator_rs::Commutator;
use lie_rs::prelude::*;
use lyndon_rs::{generators::ENotation, prelude::*};
use ordered_float::NotNan;

let e = ENotation::alphabet(2);

// Create two simple Lie series
let x = LyndonWord::try_from(vec![e[0].clone()]).unwrap();
let y = LyndonWord::try_from(vec![e[1].clone()]).unwrap();

// Coefficient types need `Eq + Hash + Ord`; wrap raw floats in `NotNan`.
let one = NotNan::new(1.0).unwrap();
let x_series = LieSeries::<ENotation, NotNan<f64>>::new(vec![x.clone()], vec![one.clone()]);
let y_series = LieSeries::<ENotation, NotNan<f64>>::new(vec![y.clone()], vec![one]);

// Compute the commutator [x, y]
let commutator_series = x_series.commutator(&y_series);
println!("Commutator [x, y]: {:?}", commutator_series);
```

## Key Types

- `LieSeriesGenerator`: Trait for generating Lie series from various inputs
- `BchSeriesGenerator`: Specialized generator for Baker-Campbell-Hausdorff computations
- `LieSeries<T, U>`: Generic container for Lie algebraic series
- `RootedTree<T>`: Trees used for BCH computations
