# signature-rs

A Rust library for computing path signatures and log signatures from multidimensional data.

## Overview

This crate implements signature methods for analyzing paths and time series data:

- **Log Signatures**: Computation of logarithmic signatures using Baker-Campbell-Hausdorff theory
- **Path Analysis**: Tools for extracting geometric and topological features from paths
- **Multidimensional Support**: Handle paths in arbitrary dimensions via `ndarray` integration
- **Configurable Precision**: Support for various numeric types including rationals and floats

## Key Features

### Log Signature Computation
- `LogSignatureBuilder<T>`: Configurable builder for log signature computation
- Degree-based truncation for computational efficiency
- Integration with Lyndon bases for coefficient ordering

### Data Integration
- `ndarray` support for efficient multidimensional array operations
- Generic numeric types supporting both exact and approximate arithmetic
- Flexible input formats for path data

## Dependencies

- `ndarray`: Multidimensional array operations
- `num-traits`: Generic numeric traits
- `ordered-float`: Floating-point arithmetic
- `lie-rs`: Lie algebra computations
- `lyndon-rs`: Lyndon word foundations
- `commutator-rs`: Commutator operations

## Development Dependencies

- `rstest`: Parametric testing framework
- `num-rational`: Exact rational arithmetic
