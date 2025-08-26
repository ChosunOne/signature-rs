# lyndon-py

Python bindings for Lyndon word operations and canonical bases for free Lie algebras.

## Installation

```bash
pip install lyndon-py
```

## Usage

### Basic Lyndon Word Operations

```python
from lyndon_py import LyndonWord

# Create a Lyndon word
word = LyndonWord([0, 1, 2])
print(f"Word: {word}")
print(f"Length: {len(word)}")
print(f"Is empty: {word.is_empty()}")

# Check ordering relationships
word1 = LyndonWord([0, 1])
word2 = LyndonWord([0, 2])
print(f"{word1} < {word2}: {word1 < word2}")

# Lyndon word multiplication (concatenation)
combined = word1 * word2
print(f"Concatenation: {combined}")
```

### Lyndon Basis Generation

```python
from lyndon_py import LyndonBasis, Sort

# Create basis with lexicographical ordering
basis = LyndonBasis(alphabet_size=3, sort=Sort.Lexicographical)

# Generate all Lyndon words up to length 4
words = basis.generate_basis(max_length=4)
for i, word in enumerate(words):
    print(f"Basis element {i}: {word}")

# Get word count per degree
counts = basis.number_of_words_per_degree(max_degree=5)
print(f"Words per degree: {counts}")
```

### Ordering Types

```python
from lyndon_py import LyndonBasis, Sort

# Lexicographical ordering
lex_basis = LyndonBasis(2, Sort.Lexicographical)
lex_words = lex_basis.generate_basis(3)
print("Lexicographical:", [str(w) for w in lex_words])

# Topological ordering
top_basis = LyndonBasis(2, Sort.Topological) 
top_words = top_basis.generate_basis(3)
print("Topological:", [str(w) for w in top_words])
```

### Factorization

```python
from lyndon_py import LyndonWord

# Standard factorization
word = LyndonWord([0, 1, 1, 2])
left_factor, right_factor = word.factorize()
print(f"Factorization of {word}: ({left_factor}, {right_factor})")

# All right factors
word = LyndonWord([0, 1, 2, 2])
right_factors = word.right_factors()
print(f"Right factors of {word}: {[str(rf) for rf in right_factors]}")
```

### Goldberg Partitions

```python
from lyndon_py import LyndonWord

# Compute Goldberg partition
word = LyndonWord([0, 1, 0, 1])
goldberg = word.goldberg()
print(f"Goldberg partition of {word}: {goldberg}")
```

## API Reference

### LyndonWord

```python
LyndonWord(letters: list[int])
```

Creates a Lyndon word from a list of integers.

**Properties:**
- `letters: list[int]`: The letter sequence

**Methods:**
- `__len__() -> int`: Length of the word
- `is_empty() -> bool`: Check if word is empty  
- `__mul__(other: LyndonWord) -> LyndonWord`: Concatenation
- `factorize() -> tuple[LyndonWord, LyndonWord]`: Standard factorization
- `right_factors() -> list[LyndonWord]`: All proper right factors
- `goldberg() -> list[int]`: Goldberg partition

**Comparison operators:** `<`, `<=`, `>`, `>=`, `==`, `!=`

### LyndonBasis

```python
LyndonBasis(alphabet_size: int, sort: Sort)
```

Factory for generating Lyndon word bases.

**Properties:**
- `alphabet_size: int`: Size of the alphabet
- `sort: Sort`: Ordering scheme

**Methods:**
- `generate_basis(max_length: int) -> list[LyndonWord]`: Generate words up to given length
- `number_of_words_per_degree(max_degree: int) -> list[int]`: Count words at each degree

### Sort

```python
class Sort(Enum):
    Lexicographical  # Dictionary ordering
    Topological      # Degree-lexicographical ordering
```