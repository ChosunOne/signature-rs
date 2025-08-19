# lie-rs

A Rust library for Lie algebra computations and Baker-Campbell-Hausdorff series.

## Overview

This crate implements algebraic structures and computations including:

- **Baker-Campbell-Hausdorff (BCH) Series**: Computation of the BCH formula for non-commuting elements
    - Adapted from https://github.com/HaraldHofstaetter/BCH
- **Lie Series**: Representation and manipulation of formal power series in Lie algebras

## Key Components

### Series Generation
- `LieSeriesGenerator`: Trait for generating Lie series from various inputs
- `BchSeriesGenerator`: Specialized generator for Baker-Campbell-Hausdorff computations
- `LieSeries<T, U>`: Generic container for Lie algebraic series
- Support for rational, floating-point, and integer arithmetic

## Features

- **Default**: Core functionality
- **Progress** (optional): Progress bar support via `indicatif`

## Dependencies

- `commutator-rs`: Commutator operations
- `lyndon-rs`: Lyndon word foundations
- `num-traits`: Generic numeric operations
- `rayon`: Parallel computation support
- `indicatif`: Progress indicators (optional)

## Development Dependencies

- `rstest`: Parametric testing
- `num-rational`: Rational number arithmetic
- `ordered-float`: Floating-point comparisons
