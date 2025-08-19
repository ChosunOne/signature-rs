# signature-rs-bin

A command-line tool for computing log signatures from path data in various formats.

## Overview

This binary provides a practical interface to the signature-rs library ecosystem, enabling users to:

- **Load Path Data**: Read multidimensional path data from CSV, JSON, JSONL, and Parquet files
- **Compute Log Signatures**: Calculate log signatures up to a specified degree using Baker-Campbell-Hausdorff theory
- **Export Results**: Output coefficients in multiple formats (CSV, JSON, JSONL, Parquet)
- **Flexible Precision**: Support for both 32-bit and 64-bit floating-point computation

## Usage

```bash
signature-rs-bin [OPTIONS] -p <PATH>
```

### Arguments

- `-d <NUM_DIMENSIONS>`: Number of dimensions in the path (auto-detected if not specified)
- `-k <MAX_DEGREE>`: Maximum degree for log signature computation (default: 3)
- `-t <DATA_TYPE>`: Coefficient data type - `f32` or `f64` (default: f32)
- `-p <PATH>`: Input file path containing the path data
- `-o <OUTPUT>`: Output file path (prints to stdout if not specified)
- `-f <OUTPUT_TYPE>`: Output format - `csv`, `json`, `jsonl`, or `parquet` (default: csv)

### Example

```bash
# Compute log signature from CSV path data, output to JSON
signature-rs-bin -p path_data.csv -o results.json -f json -k 5

# Use 64-bit precision with Parquet output
signature-rs-bin -p data.parquet -t f64 -o output.parquet -f parquet
```

## Supported Input Formats

- **CSV**: Comma-separated values (no header expected)
- **JSON**: Standard JSON format
- **JSONL**: JSON Lines format
- **Parquet**: Apache Parquet columnar format
