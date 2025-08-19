# commutator-rs

A Rust library for commutator operations and formal expressions in algebraic computation.

## Overview

This crate implements:

- **Commutator Operations**: The fundamental `[A, B] = AB - BA` operation
- **Commutator Terms**: Structured representation of algebraic expressions involving commutators
- **Formal Indeterminates**: Abstract algebraic elements for symbolic computation
- **Expression Trees**: Hierarchical representation of nested commutator expressions

## Key Features

- Generic commutator trait implementation for types supporting multiplication and subtraction
- Macro support for convenient commutator notation: `comm!(a, b)`
- Coefficient-based algebraic expressions
- Integration with Lyndon word ordering from `lyndon-rs`

## Core Types

- `Commutator`: Trait defining the commutator operation
- `CommutatorTerm<T, U>`: Represents either atomic terms or nested commutator expressions
- Support for both owned and borrowed operands

## Dependencies

- `lyndon-rs`: Lyndon word foundations
- `ordered-float`: Floating-point arithmetic
- `num-traits`: Numeric trait abstractions
- `rstest`: Testing framework (dev dependency)
