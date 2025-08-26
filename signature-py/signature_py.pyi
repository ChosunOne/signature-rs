"""Path signature computation for rough path theory and machine learning."""

from typing import Self
from lie_py import LieSeries
from numpy import float32
from numpy.typing import NDArray

class LogSignature:
    """Computed log signature of a multidimensional path."""
    
    def __getitem__(self, idx: int) -> float:
        """Get the signature coefficient at a specific index."""
        ...
    
    def __setitem__(self, idx: int, coefficient: float) -> None:
        """Set the signature coefficient at a specific index."""
        ...
    
    @property
    def series(self) -> LieSeries:
        """Get the underlying Lie series representation."""
        ...
    
    @property
    def bch_series(self) -> LieSeries:
        """Get the Baker-Campbell-Hausdorff series used for concatenation."""
        ...
    
    def concatenate(self, other: Self) -> Self:
        """Concatenate two log signatures using the Baker-Campbell-Hausdorff formula."""
        ...
    
    def concatenate_assign(self, other: Self) -> None:
        """In-place concatenation of log signatures."""
        ...

class LogSignatureBuilder:
    """Factory class for building log signatures with configurable parameters."""
    
    def __init__(
        self, max_degree: int | None = None, num_dimensions: int | None = None
    ) -> None:
        """Create a new log signature builder."""
        ...
    
    @property
    def max_degree(self) -> int:
        """Get the maximum degree of signature computation."""
        ...
    
    @max_degree.setter
    def max_degree(self, max_degree: int) -> None:
        """Set the maximum degree of signature computation."""
        ...
    
    @property
    def num_dimensions(self) -> int:
        """Get the number of path dimensions."""
        ...
    
    @num_dimensions.setter
    def num_dimensions(self, num_dimensions: int) -> None:
        """Set the number of path dimensions."""
        ...
    
    def build(self) -> LogSignature:
        """Create an empty log signature with current builder settings."""
        ...
    
    def build_from_path(self, path: NDArray[float32]) -> LogSignature:
        """Compute the log signature of a multidimensional path."""
        ...