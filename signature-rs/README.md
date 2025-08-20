# signature-rs

A Rust library for computing log signatures from path data.

## Overview

- Log signature computation from paths
- Multidimensional path support via `ndarray`
- Generic numeric types (rationals, floats, integers)

## Quick Start

Add this to your `Cargo.toml`:

```toml
[dependencies]
signature-rs = "0.1.0"
ndarray = "0.16"
```

### Basic Example

```rust
use signature_rs::prelude::*;
use ndarray::array;
use ordered_float::NotNan;

// Create a path in 2D space
let path = array![[0.0, 0.0], [1.0, 0.5], [2.0, 1.0], [3.0, 0.0]];
let path = path.mapv(|v| NotNan::new(v).expect("value to be a number"));

// Build and compute log signature up to degree 3
let builder = LogSignatureBuilder::<ENotation>::new()
    .with_num_dimensions(2)
    .with_max_degree(3);

let log_sig = builder.build_from_path(&path.view());

// Access the computed log signature
println!("Log signature computed with {} terms", log_sig.series.coefficients.len());
```

### Advanced Example with Custom Types

```rust
use signature_rs::prelude::*;
use ndarray::Array2;
use num_rational::Ratio;

// Use exact rational arithmetic for precise computation
type Rational = Ratio<i64>;

// Create a path with rational coordinates
let mut path = Array2::<Rational>::zeros((4, 2));
path[[0, 0]] = Ratio::new(0, 1);
path[[0, 1]] = Ratio::new(0, 1);
path[[1, 0]] = Ratio::new(1, 1);
path[[1, 1]] = Ratio::new(1, 2);
path[[2, 0]] = Ratio::new(2, 1);
path[[2, 1]] = Ratio::new(1, 1);

// Compute log signature with exact arithmetic
let builder = LogSignatureBuilder::<ENotation>::new()
    .with_num_dimensions(2)
    .with_max_degree(4);

let log_sig = builder.build_from_path(&path.view());

// The result contains exact rational coefficients
for (word, coeff) in log_sig.series.basis.iter().zip(log_sig.series.coefficients.iter()) {
    println!("Word: {:?}, Coefficient: {}", word, coeff);
}
```


## Key Types

- `LogSignatureBuilder<T>`: Configurable builder for log signature computation
- `LogSignature<T, U>`: Contains computed log signature data

