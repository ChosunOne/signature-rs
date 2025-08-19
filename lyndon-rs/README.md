# lyndon-rs

A Rust library for working with Lyndon words and bases in combinatorial mathematics.

## Overview

This crate provides functionality for:

- **Lyndon Words**: Working with lexicographically minimal strings in their conjugacy class
- **Lyndon Basis**: Generating and managing ordered bases using Lyndon words
- **Generator Support**: Abstract generator types for constructing Lyndon words
- **Sorting Methods**: Both lexicographical and topological sorting of Lyndon words

## Key Types

- `LyndonWord<T>`: Represents a Lyndon word over alphabet type T
- `LyndonBasis<T>`: Manages a basis of Lyndon words with configurable sorting
- `Generator`: Trait for types that can be used as alphabet elements

## Sorting Options

- **Lexicographical**: Standard dictionary ordering
- **Topological**: Order based on structural properties

## Dependencies

- `thiserror`: Error handling
- `rstest`: Testing framework (dev dependency)

## Usage

This crate is primarily used as a foundation for other algebraic computation libraries in the signature-rs workspace, particularly for constructing bases in free Lie algebras and related mathematical structures.