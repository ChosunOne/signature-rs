# lyndon-rs

A Rust library for working with Lyndon words.

## Overview

- Lyndon word generation and validation
- Lyndon basis construction
- Generator types for alphabet elements
- Multiple sorting methods

## Quick Start

Add this to your `Cargo.toml`:

```toml
[dependencies]
lyndon-rs = "0.1.0"
```

### Basic Example

```rust
use lyndon_rs::prelude::*;

// Create a Lyndon basis with 3 generators up to degree 4
let basis = LyndonBasis::<ENotation>::new(3, Sort::Lexicographical);

// Generate all Lyndon words up to degree 4
let words = basis.generate_basis(4);

println!("Generated {} Lyndon words", words.len());
for word in &words[..5] {  // Show first 5
    println!("Word: {:?}, Length: {}", word, word.len());
}
```

### Working with Custom Generators

```rust
use lyndon_rs::prelude::*;

// Define a custom generator type
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct Letter(char);

impl Generator for Letter {
    fn generate_set(size: usize) -> Vec<Self> {
        (0..size).map(|i| Letter((b'a' + i as u8) as char)).collect()
    }
}

// Create basis with custom generators
let basis = LyndonBasis::<Letter>::new(4, Sort::Lexicographical);  // a, b, c, d

let words = basis.generate_basis(3);
for word in words {
    println!("{:?}", word);
}
```

### Factorizing Lyndon Words

```rust
use lyndon_rs::prelude::*;

// Create a specific Lyndon word
let word = LyndonWord::try_from(vec![ENotation::new(1), ENotation::new(2), ENotation::new(1)])
    .expect("Valid Lyndon word");

// Factor it into right factors
let factors = word.right_factors();
println!("Right factors of {:?}:", word);
for factor in factors {
    println!("  {:?}", factor);
}

// Get Goldberg representation for algebraic operations
let goldberg_rep = word.goldberg();
println!("Goldberg representation: {:?}", goldberg_rep);

// Factor into canonical factorization (v, w)
if word.len() > 1 {
    let (left, right) = word.factorize();
    println!("Factorization: ({:?}, {:?})", left, right);
}
```

### Constructing Words and Checking Properties

```rust
use lyndon_rs::prelude::*;

// Check if a word is Lyndon
let generators = vec![ENotation::new(1), ENotation::new(2), ENotation::new(1), ENotation::new(3)];
match LyndonWord::try_from(generators.clone()) {
    Ok(lyndon_word) => {
        println!("{:?} is a valid Lyndon word", lyndon_word);
        println!("Length: {}", lyndon_word.len());
        println!("Degree: {}", lyndon_word.degree());
    }
    Err(e) => {
        println!("{:?} is not a Lyndon word: {}", generators, e);
    }
}

// Concatenate two Lyndon words (may not result in a valid Lyndon word)
let word1 = LyndonWord::try_from(vec![ENotation::new(1)]).unwrap();
let word2 = LyndonWord::try_from(vec![ENotation::new(2), ENotation::new(3)]).unwrap();

match word1 * word2 {
    Ok(concatenated) => println!("Concatenated word: {:?}", concatenated),
    Err(e) => println!("Concatenation failed: {}", e),
}
```


## Key Types

- `LyndonWord<T>`: Represents a Lyndon word over alphabet type T
- `LyndonBasis<T>`: Manages a basis of Lyndon words with configurable sorting
- `Generator`: Trait for types that can be used as alphabet elements
- `ENotation`: Standard e-notation generators (e1, e2, e3, ...)
- `Sort`: Sorting methods (`Lexicographical`, `Topological`)