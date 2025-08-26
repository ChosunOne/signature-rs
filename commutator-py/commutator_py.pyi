"""Commutator operations and formal indeterminates for Lie algebra computations."""

from typing import Self, override
from lyndon_py import LyndonWord

class FormalIndeterminate:
    """Symbolic representation of formal indeterminate expressions."""
    
    def __init__(self, coefficient: float, symbols: list[int]) -> None:
        """Create a new formal indeterminate."""
        ...
    
    @classmethod
    def from_commutator(cls, commutator: CommutatorTerm) -> list[Self]:
        """Convert a commutator term to its formal indeterminate representation."""
        ...
    
    def __mul__(self, other: float | Self) -> Self:
        """Multiply by a scalar or another indeterminate."""
        ...
    
    def __rmul__(self, other: float | Self) -> Self:
        """Right multiplication by scalar or indeterminate."""
        ...
    
    def __neg__(self) -> Self:
        """Negate the formal indeterminate."""
        ...
    
    @override
    def __eq__(self, other: object) -> bool:
        """Test equality with another formal indeterminate."""
        ...
    
    @override
    def __hash__(self) -> int:
        """Return hash value for use in sets and dictionaries."""
        ...
    
    @property
    def coefficient(self) -> float:
        """Get the numerical coefficient."""
        ...
    
    @coefficient.setter
    def coefficient(self, coefficient: float) -> None:
        """Set the numerical coefficient."""
        ...
    
    @property
    def symbols(self) -> list[int]:
        """Get the generator symbols."""
        ...
    
    @symbols.setter
    def symbols(self, symbols: list[int]) -> None:
        """Set the generator symbols."""
        ...

class CommutatorTerm:
    """Commutator expression with nested bracket structure."""
    
    def __init__(
        self,
        coefficient: float,
        atom: int | None = None,
        left: Self | None = None,
        right: Self | None = None,
    ) -> None:
        """Create a new commutator term."""
        ...
    
    def __mul__(self, other: float) -> Self:
        """Multiply the commutator term by a scalar."""
        ...
    
    def __rmul__(self, other: float) -> Self:
        """Right multiplication by scalar."""
        ...
    
    def __neg__(self) -> Self:
        """Negate the commutator term."""
        ...
    
    def __lt__(self, other: Self) -> bool:
        """Compare commutator terms using canonical ordering."""
        ...
    
    def __le__(self, other: Self) -> bool:
        """Less than or equal comparison."""
        ...
    
    def __gt__(self, other: Self) -> bool:
        """Greater than comparison."""
        ...
    
    def __ge__(self, other: Self) -> bool:
        """Greater than or equal comparison."""
        ...
    
    @override
    def __eq__(self, other: object) -> bool:
        """Test equality with another commutator term."""
        ...
    
    @override
    def __hash__(self) -> int:
        """Return hash value for use in sets and dictionaries."""
        ...
    
    @property
    def atom(self) -> int | None:
        """Get the atomic generator index."""
        ...
    
    @atom.setter
    def atom(self, atom: int) -> None:
        """Set the atomic generator index."""
        ...
    
    def degree(self) -> int:
        """Get the degree (depth) of the commutator expression."""
        ...
    
    @property
    def coefficient(self) -> float:
        """Get the numerical coefficient."""
        ...
    
    @coefficient.setter
    def coefficient(self, coefficient: float) -> None:
        """Set the numerical coefficient."""
        ...
    
    def commutator(self, rhs: Self) -> Self:
        """Compute the commutator [self, rhs] = self·rhs - rhs·self."""
        ...
    
    @classmethod
    def from_lyndon_word(cls, lyndon_word: LyndonWord) -> Self:
        """Create a commutator term from a Lyndon word."""
        ...
    
    def formal_indeterminates(self) -> list[FormalIndeterminate]:
        """Convert to formal indeterminate representation."""
        ...
    
    def is_zero(self) -> bool:
        """Check if the commutator term represents zero."""
        ...
    
    def jacobi_identity(self) -> tuple[Self, Self] | None:
        """Compute terms for Jacobi identity verification."""
        ...
    
    @property
    def left(self) -> Self | None:
        """Get the left operand for composite terms."""
        ...
    
    @left.setter
    def left(self, left: Self) -> None:
        """Set the left operand."""
        ...
    
    def lyndon_basis_decomposition(self, lyndon_basis_set: set[Self]) -> list[Self]:
        """Decompose the commutator term using a Lyndon basis."""
        ...
    
    def lyndon_sort(self) -> None:
        """Sort the commutator term into Lyndon canonical form."""
        ...
    
    @property
    def right(self) -> Self | None:
        """Get the right operand for composite terms."""
        ...
    
    @right.setter  
    def right(self, right: Self) -> None:
        """Set the right operand."""
        ...
    
    def unit(self) -> Self:
        """Create a unit commutator term with the same structure."""
        ...