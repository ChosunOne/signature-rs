# signature-py

Python bindings for path signature computation.

## Installation

```bash
pip install signature-py
```

## Usage

### Basic Log Signature Computation

```python
import numpy as np
from signature_py import LogSignatureBuilder

# Create a simple 2D path
path = np.array([
    [0.0, 0.0],
    [1.0, 0.5],
    [2.0, 1.0],
    [1.5, 2.0]
], dtype=np.float32)

# Build log signature with degree 3
builder = LogSignatureBuilder(max_degree=3, num_dimensions=2)
log_signature = builder.build_from_path(path)

# Access signature coefficients
coefficients = log_signature.coefficients
print(f"Log signature coefficients: {coefficients}")

# Get individual coefficient values
for i in range(len(coefficients)):
    print(f"Coefficient {i}: {log_signature[i]}")
```

### Path Concatenation

```python
import numpy as np
from signature_py import LogSignatureBuilder

# Create two path segments
path1 = np.array([[0., 0.], [1., 1.]], dtype=np.float32)
path2 = np.array([[1., 1.], [2., 0.]], dtype=np.float32) 

builder = LogSignatureBuilder(max_degree=4, num_dimensions=2)

# Compute signatures separately
sig1 = builder.build_from_path(path1)
sig2 = builder.build_from_path(path2)

# Concatenate using Baker-Campbell-Hausdorff formula
concatenated_sig = sig1.concatenate(sig2)
print(f"Concatenated signature: {concatenated_sig.coefficients}")

# In-place concatenation
sig1.concatenate_assign(sig2)
print(f"In-place result: {sig1.coefficients}")
```

### Working with High-Dimensional Data

```python
import numpy as np
from signature_py import LogSignatureBuilder

# Generate 5D time series data
np.random.seed(42)
n_points, n_dims = 100, 5
path_data = np.cumsum(np.random.randn(n_points, n_dims), axis=0).astype(np.float32)

# Compute signature with degree 3
builder = LogSignatureBuilder(max_degree=3, num_dimensions=5)
signature = builder.build_from_path(path_data)

print(f"Signature dimension: {len(signature.coefficients)}")
print(f"First 10 coefficients: {signature.coefficients[:10]}")
```

### Empty Signature Construction

```python
from signature_py import LogSignatureBuilder

# Create empty signature for initialization
builder = LogSignatureBuilder(max_degree=4, num_dimensions=3)
empty_sig = builder.build()

print(f"Empty signature length: {len(empty_sig.coefficients)}")
print(f"All zeros: {all(c == 0.0 for c in empty_sig.coefficients)}")
```

## API Reference

### LogSignature

#### Properties

- `coefficients: list[float]`: The underlying Lie series coefficients
- `series: LieSeries`: The underlying Lie series representation
- `bch_series: LieSeries`: BCH series used for concatenation operations

#### Methods

- `__getitem__(idx: int) -> float`: Get coefficient at specific index
- `__setitem__(idx: int, coefficient: float) -> None`: Set coefficient at specific index
- `concatenate(other: LogSignature) -> LogSignature`: Concatenate two signatures using BCH
- `concatenate_assign(other: LogSignature) -> None`: In-place concatenation

### LogSignatureBuilder

#### Constructor

```python
LogSignatureBuilder(max_degree: int | None = None, num_dimensions: int | None = None)
```

Creates a signature builder with specified parameters.

#### Properties

- `max_degree: int`: Maximum degree of signature terms
- `num_dimensions: int`: Number of dimensions in path data

#### Methods

- `build() -> LogSignature`: Create empty signature with current settings
- `build_from_path(path: NDArray[float32]) -> LogSignature`: Compute signature from path data

**Path Requirements:**
- `path` must be a numpy array with shape `(n_points, n_dimensions)`
- `path` must have dtype `np.float32`
- `path` must have at least 2 points
- `n_dimensions` must match the builder's `num_dimensions`
