# commutator-py

Python bindings for commutator operations and formal indeterminates for Lie algebra computations.

## Installation

```bash
pip install commutator-py
```

## Usage

### Basic Commutator Operations

```python
from commutator_py import CommutatorTerm

# Create atomic terms
a = CommutatorTerm(1.0, atom=0)
b = CommutatorTerm(1.0, atom=1)

# Compute commutator [a, b]
commutator_ab = a.commutator(b)
print(f"[a, b] = {commutator_ab}")

# Create nested commutators: [[a, b], c]
c = CommutatorTerm(1.0, atom=2)
nested = commutator_ab.commutator(c)
print(f"[[a, b], c] = {nested}")
```

### Working with Lyndon Words

```python
from commutator_py import CommutatorTerm
from lyndon_py import LyndonWord

# Convert Lyndon word to commutator term
lyndon_word = LyndonWord([0, 1, 1])
commutator = CommutatorTerm.from_lyndon_word(lyndon_word)
print(f"Commutator from Lyndon word: {commutator}")
```

### Formal Indeterminates

```python
from commutator_py import FormalIndeterminate, CommutatorTerm

# Create a commutator term
term = CommutatorTerm(2.0, left=CommutatorTerm(1.0, 0), right=CommutatorTerm(1.0, 1))

# Convert to formal indeterminates
indeterminates = FormalIndeterminate.from_commutator(term)
for ind in indeterminates:
    print(f"Coefficient: {ind.coefficient}, Symbols: {ind.symbols}")
```

## API Reference

### CommutatorTerm

#### Constructors

```python
# Atomic term
CommutatorTerm(coefficient: float, atom: int)

# Composite term
CommutatorTerm(coefficient: float, left: CommutatorTerm, right: CommutatorTerm)

# From Lyndon word
CommutatorTerm.from_lyndon_word(lyndon_word: LyndonWord) -> CommutatorTerm
```

#### Methods

- `commutator(other: CommutatorTerm) -> CommutatorTerm`: Compute `[self, other]`
- `degree() -> int`: Get the degree (depth) of the commutator
- `is_zero() -> bool`: Check if the term represents zero
- `formal_indeterminates() -> list[FormalIndeterminate]`: Convert to formal representation
- `lyndon_basis_decomposition(basis: set[CommutatorTerm]) -> list[CommutatorTerm]`: Decompose using Lyndon basis

#### Properties

- `coefficient: float`: Numerical coefficient
- `atom: int | None`: Atomic generator index (for leaf terms)
- `left: CommutatorTerm | None`: Left operand (for composite terms)
- `right: CommutatorTerm | None`: Right operand (for composite terms)

### FormalIndeterminate

#### Constructor

```python
FormalIndeterminate(coefficient: float, symbols: list[int])
```

#### Methods

- `from_commutator(commutator: CommutatorTerm) -> list[FormalIndeterminate]`: Convert commutator to formal representation

#### Properties

- `coefficient: float`: Numerical coefficient
- `symbols: list[int]`: Ordered list of generator indices