"""Lyndon word operations and canonical bases for free Lie algebras."""

from enum import Enum, auto
from typing import Self, override

class Sort(Enum):
    """Ordering schemes for Lyndon word generation."""
    Lexicographical = auto()
    Topological = auto()

class LyndonWord:
    """A Lyndon word with validation and operations."""
    
    def __init__(self, letters: list[int]) -> None:
        """Create a Lyndon word from a list of integers."""
        ...
    
    @property
    def letters(self) -> list[int]:
        """Get the underlying letter sequence."""
        ...
    
    @letters.setter
    def letters(self, letters: list[int]) -> None:
        """Set the letter sequence."""
        ...
    
    def __len__(self) -> int:
        """Get the length of the word."""
        ...
    
    def is_empty(self) -> bool:
        """Check if the word is empty."""
        ...
    
    def __mul__(self, other: Self) -> Self:
        """Concatenate two Lyndon words."""
        ...
    
    def __rmul__(self, other: Self) -> Self:
        """Right multiplication."""
        ...
    
    def __lt__(self, other: Self) -> bool:
        """Less than comparison."""
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
        """Equality comparison."""
        ...
    
    @override
    def __ne__(self, other: object) -> bool:
        """Inequality comparison."""
        ...
    
    @override
    def __hash__(self) -> int:
        """Hash value for use in sets and dictionaries."""
        ...
    
    def factorize(self) -> tuple[Self, Self]:
        """Return the standard factorization (u, v) where word = uv."""
        ...
    
    def right_factors(self) -> list[Self]:
        """Get all proper right factors of the word."""
        ...
    
    def goldberg(self) -> list[int]:
        """Compute the Goldberg partition."""
        ...

class LyndonBasis:
    """Factory for generating Lyndon word bases."""
    
    def __init__(self, alphabet_size: int, sort: Sort) -> None:
        """Create a Lyndon basis generator."""
        ...
    
    @property
    def alphabet_size(self) -> int:
        """Get the alphabet size."""
        ...
    
    @alphabet_size.setter
    def alphabet_size(self, alphabet_size: int) -> None:
        """Set the alphabet size."""
        ...
    
    @property
    def sort(self) -> Sort:
        """Get the sort order."""
        ...
    
    @sort.setter
    def sort(self, sort: Sort) -> None:
        """Set the sort order."""
        ...
    
    def generate_basis(self, max_length: int) -> list[LyndonWord]:
        """Generate all Lyndon words up to the specified length."""
        ...
    
    def number_of_words_per_degree(self, max_degree: int) -> list[int]:
        """Get the count of words at each degree level."""
        ...