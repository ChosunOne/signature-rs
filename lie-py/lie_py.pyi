"""Lie series and Baker-Campbell-Hausdorff computations."""

from typing import Self
from commutator_py import CommutatorTerm
from lyndon_py import LyndonBasis, LyndonWord

class LieSeries:
    """Formal Lie series with algebraic operations over a Lyndon basis."""
    
    def __init__(self, basis: list[LyndonWord], coefficients: list[float]) -> None:
        """Create a new Lie series from basis elements and coefficients."""
        ...
    
    def __len__(self) -> int:
        """Get the number of basis elements in the series."""
        ...
    
    def __getitem__(self, idx: int) -> float:
        """Get the coefficient at a specific index."""
        ...
    
    def __setitem__(self, idx: int, coefficient: float) -> None:
        """Set the coefficient at a specific index."""
        ...
    
    def __add__(self, other: Self) -> Self:
        """Add two Lie series."""
        ...
    
    def __radd__(self, other: Self) -> Self:
        """Right addition."""
        ...
    
    def __iadd__(self, other: Self) -> None:
        """In-place addition."""
        ...
    
    def __sub__(self, other: Self) -> Self:
        """Subtract two Lie series."""
        ...
    
    def __rsub__(self, other: Self) -> Self:
        """Right subtraction."""
        ...
    
    def __isub__(self, other: Self) -> None:
        """In-place subtraction."""
        ...
    
    def __mul__(self, other: float) -> Self:
        """Multiply Lie series by a scalar."""
        ...
    
    def __rmul__(self, other: float) -> Self:
        """Right scalar multiplication."""
        ...
    
    def __imul__(self, other: float) -> None:
        """In-place scalar multiplication."""
        ...
    
    @property
    def basis(self) -> list[LyndonWord]:
        """Get the Lyndon word basis."""
        ...
    
    @property
    def coefficients(self) -> list[float]:
        """Get the coefficient vector."""
        ...
    
    @property
    def commutator_basis(self) -> list[CommutatorTerm]:
        """Get the commutator term representation of the basis."""
        ...
    
    def commutator(self, other: Self) -> Self:
        """Compute the Lie bracket [self, other] of two Lie series."""
        ...
    
    @property
    def max_degree(self) -> int:
        """Get the maximum degree of terms in the series."""
        ...

class BCHSeriesGenerator:
    """Generator for Baker-Campbell-Hausdorff series with exact coefficient computation."""
    
    def __init__(self, basis: LyndonBasis, max_degree: int) -> None:
        """Create a new BCH series generator."""
        ...
    
    @property
    def alphabet_size(self) -> int:
        """Get the size of the generating alphabet."""
        ...
    
    @property
    def basis(self) -> list[LyndonWord]:
        """Get the Lyndon basis used for the BCH series."""
        ...
    
    def generate_goldberg_coefficient_numerators(self) -> list[int]:
        """Generate exact Goldberg coefficient numerators for the BCH series."""
        ...
    
    def generate_lie_series(self) -> LieSeries:
        """Generate the universal Baker-Campbell-Hausdorff Lie series."""
        ...
    
    @property
    def index_of_degree(self) -> list[int]:
        """Get starting indices for each degree in the basis."""
        ...
    
    @property
    def max_degree(self) -> int:
        """Get the maximum degree used in the BCH series."""
        ...
    
    @property
    def multi_degree(self) -> list[int]:
        """Get multi-degree information for basis elements."""
        ...
    
    @property
    def left_factor(self) -> list[int]:
        """Get left factor indices for BCH factorization."""
        ...
    
    @property
    def right_factor(self) -> list[int]:
        """Get right factor indices for BCH factorization."""
        ...
    
    @property
    def word_lengths(self) -> list[int]:
        """Get the length of each basis word."""
        ...