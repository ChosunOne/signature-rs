# lie-py

Python bindings for Lie series and Baker-Campbell-Hausdorff computations.

## Installation

```bash
pip install lie-py
```

## Usage

### Basic Lie Series Operations

```python
from lie_py import LieSeries
from lyndon_py import LyndonBasis, Sort

# Create a Lyndon basis
basis = LyndonBasis(alphabet_size=2, sort=Sort.Lexicographical)
lyndon_words = basis.generate_basis(max_length=3)

# Create Lie series
x = LieSeries(lyndon_words, [1.0, 0.0, 0.5, 0.0, 0.0])
y = LieSeries(lyndon_words, [0.0, 1.0, 0.0, 0.25, 0.0])

# Compute commutator [X, Y]
commutator_xy = x.commutator(y)
print(f"Coefficients of [X,Y]: {commutator_xy.coefficients}")

# Add series: X + Y
sum_xy = x + y
print(f"Coefficients of X+Y: {sum_xy.coefficients}")
```

### Baker-Campbell-Hausdorff Series Generation

```python
from lie_py import BCHSeriesGenerator
from lyndon_py import LyndonBasis, Sort

# Set up BCH computation
basis = LyndonBasis(alphabet_size=2, sort=Sort.Topological)
max_degree = 4

# Generate BCH series structure
bch_generator = BCHSeriesGenerator(basis, max_degree)

# Get the universal BCH series
bch_series = bch_generator.generate_lie_series()
print(f"BCH series has {len(bch_series)} terms up to degree {max_degree}")

# Access coefficient numerators (Goldberg coefficients)
numerators = bch_generator.generate_goldberg_coefficient_numerators()
print(f"Goldberg numerators: {numerators}")
```

### Working with Commutator Representations

```python
from lie_py import LieSeries
from commutator_py import CommutatorTerm
from lyndon_py import LyndonBasis, Sort

# Create Lie series
basis = LyndonBasis(2, Sort.Lexicographical)
words = basis.generate_basis(3)
series = LieSeries(words, [1.0, 2.0, 0.5, -0.25, 1.5])

# Access as commutator terms
commutator_basis = series.commutator_basis
for i, term in enumerate(commutator_basis):
    coeff = series.coefficients[i]
    print(f"Term {i}: {coeff} * {term}")
```

## API Reference

### LieSeries

#### Constructor

```python
LieSeries(basis: list[LyndonWord], coefficients: list[float])
```

Creates a Lie series from basis elements and coefficients.

#### Properties

- `basis: list[LyndonWord]`: The Lyndon word basis
- `coefficients: list[float]`: Coefficient vector
- `commutator_basis: list[CommutatorTerm]`: Commutator term representation
- `max_degree: int`: Maximum degree of terms

#### Methods

- `__len__() -> int`: Number of basis elements
- `__getitem__(idx: int) -> float`: Get coefficient at index
- `__setitem__(idx: int, coefficient: float) -> None`: Set coefficient at index
- `commutator(other: LieSeries) -> LieSeries`: Compute Lie bracket `[self, other]`

#### Algebraic Operations

- `__add__(other: LieSeries) -> LieSeries`: Addition
- `__iadd__(other: LieSeries) -> None`: In-place addition
- `__sub__(other: LieSeries) -> LieSeries`: Subtraction
- `__isub__(other: LieSeries) -> None`: In-place subtraction
- `__mul__(other: float) -> LieSeries`: Scalar multiplication
- `__imul__(other: float) -> None`: In-place scalar multiplication

### BCHSeriesGenerator

#### Constructor

```python
BCHSeriesGenerator(basis: LyndonBasis, max_degree: int)
```

Creates a BCH series generator with the specified basis and maximum degree.

#### Properties

- `alphabet_size: int`: Size of the generating alphabet
- `basis: list[LyndonWord]`: Lyndon basis used for the BCH series
- `index_of_degree: list[int]`: Starting indices for each degree
- `max_degree: int`: Maximum degree of generated terms
- `multi_degree: list[int]`: Multi-degree information for basis elements
- `left_factor: list[int]`: Left factor indices for BCH factorization
- `right_factor: list[int]`: Right factor indices for BCH factorization
- `word_lengths: list[int]`: Length of each basis word

#### Methods

- `generate_goldberg_coefficient_numerators() -> list[int]`: Generate exact Goldberg coefficient numerators
- `generate_lie_series() -> LieSeries`: Generate the universal BCH Lie series