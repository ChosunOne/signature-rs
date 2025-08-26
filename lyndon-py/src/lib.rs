//! Python bindings for Lyndon word operations and canonical bases for free Lie algebras.
//!
//! # Classes
//! 
//! - **LyndonWord**: Individual Lyndon word manipulation with validation
//! - **LyndonBasis**: Factory for generating Lyndon word bases
//! - **Sort**: Ordering schemes (lexicographical and topological)

use std::{fmt::Display, hash::Hash, ops::Mul};

use lyndon_rs::{LyndonBasis, LyndonWord, LyndonWordError, Sort};
use pyo3::{exceptions::PyValueError, prelude::*, types::PyList};

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct LyndonWordPyError(LyndonWordError);

impl From<LyndonWordError> for LyndonWordPyError {
    fn from(value: LyndonWordError) -> Self {
        Self(value)
    }
}

impl From<LyndonWordPyError> for PyErr {
    fn from(value: LyndonWordPyError) -> Self {
        match value.0 {
            LyndonWordError::InvalidWord => PyValueError::new_err("Invalid lyndon word"),
            LyndonWordError::InvalidLetter => PyValueError::new_err("Invalid letter"),
        }
    }
}

/// Ordering schemes for Lyndon word generation and comparison.
///
/// Different orderings produce different Lyndon basis structures, affecting
/// the canonical representation of Lie algebra elements.
#[pyclass(name = "Sort", eq)]
#[repr(u8)]
#[derive(Copy, Clone, PartialEq, Eq)]
pub enum SortPy {
    /// Standard lexicographical (dictionary) ordering.
    /// 
    /// Words are ordered by comparing letters from left to right,
    /// using standard dictionary ordering rules.
    Lexicographical = 0,
    
    /// Topological ordering (degree-lexicographical).
    /// 
    /// Words are first ordered by length (degree), then lexicographically
    /// within each degree. This is often preferred for mathematical applications.
    Topological = 1,
}

#[pymethods]
impl SortPy {
    #[getter]
    #[must_use]
    pub fn get_value(&self) -> u8 {
        *self as u8
    }
}

impl From<SortPy> for Sort {
    fn from(value: SortPy) -> Self {
        match value {
            SortPy::Lexicographical => Self::Lexicographical,
            SortPy::Topological => Self::Topological,
        }
    }
}

impl From<Sort> for SortPy {
    fn from(value: Sort) -> Self {
        match value {
            Sort::Lexicographical => Self::Lexicographical,
            Sort::Topological => Self::Topological,
        }
    }
}

/// A Lyndon word - a string that is strictly smaller than all its proper rotations.
///
/// Lyndon words form canonical bases for free Lie algebras and are fundamental
/// in algebraic combinatorics and rough path theory. Each Lyndon word corresponds
/// to a unique basis element in the associated free Lie algebra.
///
/// # Mathematical Properties
///
/// - Every primitive string has a unique Lyndon word representative
/// - Lyndon words admit canonical factorization: w = uv where u,v are Lyndon and u < v  
/// - The set of Lyndon words provides a linear basis for the free Lie algebra
///
/// # Examples
///
/// ```python
/// from lyndon_py import LyndonWord
/// 
/// # Create a Lyndon word
/// word = LyndonWord([0, 1])
/// assert len(word) == 2
/// assert not word.is_empty()
/// ```
#[pyclass(name = "LyndonWord", eq, ord, str, frozen, hash)]
#[derive(Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct LyndonWordPy {
    pub inner: LyndonWord<u8>,
}

impl Display for LyndonWordPy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.inner.fmt(f)
    }
}

impl LyndonWordPy {
    #[must_use]
    pub fn inner(&self) -> &LyndonWord<u8> {
        &self.inner
    }
}

#[pymethods]
impl LyndonWordPy {
    /// Create a new Lyndon word from a sequence of letters.
    ///
    /// # Arguments
    ///
    /// * `letters` - A list of integers representing the letters of the word
    ///
    /// # Returns
    ///
    /// A new Lyndon word if the input sequence satisfies the Lyndon word property
    ///
    /// # Raises
    ///
    /// * `ValueError` - If the input sequence is not a valid Lyndon word
    ///
    /// # Examples
    ///
    /// ```python
    /// from lyndon_py import LyndonWord
    /// 
    /// # Valid Lyndon words
    /// word1 = LyndonWord([0])        # Single letter
    /// word2 = LyndonWord([0, 1])     # Binary Lyndon word
    /// word3 = LyndonWord([0, 1, 2])  # Ternary Lyndon word
    /// 
    /// # This would raise ValueError (not Lyndon)
    /// # word_invalid = LyndonWord([1, 0])
    /// ```
    #[new]
    #[must_use]
    pub fn new(letters: Vec<u8>) -> PyResult<Self> {
        Ok(Self {
            inner: LyndonWord::try_from(letters).or(Err(PyValueError::new_err("Invalid word")))?,
        })
    }

    #[must_use]
    pub fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }

    /// Concatenate two Lyndon words.
    ///
    /// Note: The result may not be a Lyndon word. This operation performs
    /// simple concatenation without preserving the Lyndon property.
    ///
    /// # Arguments
    ///
    /// * `other` - The Lyndon word to concatenate with this one
    ///
    /// # Returns
    ///
    /// A new Lyndon word representing the concatenation
    ///
    /// # Raises
    ///
    /// * `ValueError` - If the concatenation fails
    #[must_use]
    pub fn __mul__(&self, other: &Self) -> PyResult<Self> {
        Ok(Self {
            inner: (&self.inner)
                .mul(&other.inner)
                .map_err(LyndonWordPyError::from)?,
        })
    }

    /// Get the letters of the Lyndon word.
    ///
    /// # Returns
    ///
    /// A list of integers representing the letters in the word
    #[getter]
    #[must_use]
    pub fn get_letters(&self, py: Python<'_>) -> PyResult<Py<PyList>> {
        let result = PyList::new(py, &self.inner.letters)
            .map(Bound::unbind)?
            .clone_ref(py);
        PyResult::Ok(result)
    }

    /// Get the length of the Lyndon word.
    ///
    /// # Returns
    ///
    /// The number of letters in the word
    #[must_use]
    pub fn __len__(&self) -> usize {
        self.inner.len()
    }

    /// Check if the Lyndon word is empty.
    ///
    /// # Returns
    ///
    /// True if the word contains no letters, False otherwise
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// Compute the Goldberg partition of the Lyndon word.
    ///
    /// The Goldberg partition is used in Baker-Campbell-Hausdorff formula
    /// computations and provides structural information about the word.
    ///
    /// # Returns
    ///
    /// A list of integers representing the Goldberg partition
    ///
    /// # Examples
    ///
    /// ```python
    /// from lyndon_py import LyndonWord
    /// 
    /// word = LyndonWord([0, 1, 1])
    /// partition = word.goldberg()
    /// assert partition == [2, 1]  # Indicates the structure of the word
    /// ```
    #[must_use]
    pub fn goldberg(&self) -> Vec<usize> {
        self.inner.goldberg()
    }

    /// Get all proper right factors of the Lyndon word.
    ///
    /// Right factors are suffixes of the word that are also Lyndon words.
    /// This is useful for understanding the factorization structure.
    ///
    /// # Returns
    ///
    /// A list of Lyndon words representing all proper right factors
    ///
    /// # Examples
    ///
    /// ```python
    /// from lyndon_py import LyndonWord
    /// 
    /// word = LyndonWord([0, 1, 2, 2])
    /// factors = word.right_factors()
    /// # Returns all right factors that are valid Lyndon words
    /// ```
    #[must_use]
    pub fn right_factors(&self) -> Vec<Self> {
        self.inner
            .right_factors()
            .into_iter()
            .map(|inner| Self { inner })
            .collect()
    }

    /// Compute the standard factorization of the Lyndon word.
    ///
    /// Every Lyndon word w can be uniquely written as w = uv where
    /// u and v are Lyndon words and u < v in the chosen ordering.
    ///
    /// # Returns
    ///
    /// A tuple (u, v) where u and v are the factors of the factorization
    ///
    /// # Examples
    ///
    /// ```python
    /// from lyndon_py import LyndonWord
    /// 
    /// word = LyndonWord([0, 1, 1, 2])
    /// left, right = word.factorize()
    /// assert left < right  # Standard factorization property
    /// ```
    #[must_use]
    pub fn factorize(&self) -> (Self, Self) {
        let (a, b) = self.inner.factorize();
        (Self { inner: a }, Self { inner: b })
    }
}

/// Factory for generating Lyndon word bases with configurable ordering.
///
/// A Lyndon basis provides a systematic way to generate all Lyndon words
/// up to a given length over a specified alphabet, using either lexicographical
/// or topological ordering.
///
/// # Usage Patterns
///
/// - Generate complete Lyndon bases for free Lie algebra computations
/// - Count Lyndon words per degree for dimensional analysis
/// - Support both lexicographical and topological orderings
///
/// # Examples
///
/// ```python
/// from lyndon_py import LyndonBasis, Sort
/// 
/// # Create basis with lexicographical ordering
/// basis = LyndonBasis(alphabet_size=2, sort=Sort.Lexicographical)
/// words = basis.generate_basis(max_length=3)
/// 
/// # Count words per degree
/// counts = basis.number_of_words_per_degree(5)
/// ```
#[pyclass(name = "LyndonBasis")]
#[derive(Clone, Debug, Copy, Default)]
pub struct LyndonBasisPy {
    inner: LyndonBasis<u8>,
}

#[pymethods]
impl LyndonBasisPy {
    /// Create a new Lyndon basis factory.
    ///
    /// # Arguments
    ///
    /// * `alphabet_size` - The number of letters in the alphabet (must be > 0)
    /// * `sort` - The ordering scheme to use (lexicographical or topological)
    ///
    /// # Returns
    ///
    /// A new Lyndon basis factory configured with the specified parameters
    ///
    /// # Examples
    ///
    /// ```python
    /// from lyndon_py import LyndonBasis, Sort
    /// 
    /// # Binary alphabet with lexicographical ordering
    /// basis_lex = LyndonBasis(2, Sort.Lexicographical)
    /// 
    /// # Ternary alphabet with topological ordering
    /// basis_top = LyndonBasis(3, Sort.Topological)
    /// ```
    #[new]
    #[must_use]
    pub fn new(alphabet_size: usize, sort: SortPy) -> Self {
        Self {
            inner: LyndonBasis::new(alphabet_size, sort.into()),
        }
    }

    #[getter]
    #[must_use]
    pub fn get_alphabet_size(&self) -> usize {
        self.inner.alphabet_size
    }

    #[setter]
    pub fn set_alphabet_size(&mut self, alphabet_size: usize) {
        self.inner.alphabet_size = alphabet_size;
    }

    #[getter]
    #[must_use]
    pub fn get_sort(&self) -> SortPy {
        self.inner.sort.into()
    }

    #[setter]
    pub fn set_sort(&mut self, sort: SortPy) {
        self.inner.sort = sort.into();
    }

    /// Count the number of Lyndon words at each degree up to a maximum.
    ///
    /// This provides dimensional information about the free Lie algebra
    /// and is useful for understanding the growth of the basis.
    ///
    /// # Arguments
    ///
    /// * `max_degree` - Maximum degree to compute counts for
    ///
    /// # Returns
    ///
    /// A list where the i-th element is the number of Lyndon words of degree i+1
    ///
    /// # Examples
    ///
    /// ```python
    /// from lyndon_py import LyndonBasis, Sort
    /// 
    /// basis = LyndonBasis(2, Sort.Lexicographical)
    /// counts = basis.number_of_words_per_degree(4)
    /// # counts[0] = number of words of degree 1
    /// # counts[1] = number of words of degree 2
    /// # etc.
    /// ```
    #[must_use]
    pub fn number_of_words_per_degree(&self, max_degree: usize) -> Vec<usize> {
        self.inner.number_of_words_per_degree(max_degree)
    }

    /// Generate all Lyndon words up to a specified maximum length.
    ///
    /// This is the primary method for creating Lyndon bases. The words are
    /// generated in the order specified by the Sort parameter.
    ///
    /// # Arguments
    ///
    /// * `max_length` - Maximum length of Lyndon words to generate
    ///
    /// # Returns
    ///
    /// A list of Lyndon words ordered according to the basis configuration
    ///
    /// # Examples
    ///
    /// ```python
    /// from lyndon_py import LyndonBasis, Sort
    /// 
    /// basis = LyndonBasis(2, Sort.Topological)
    /// words = basis.generate_basis(3)
    /// 
    /// # Words will be in topological order:
    /// # [0], [1], [0,1], [0,0,1], [0,1,1]
    /// ```
    #[must_use]
    pub fn generate_basis(&self, max_length: usize) -> Vec<LyndonWordPy> {
        self.inner
            .generate_basis(max_length)
            .into_iter()
            .map(|inner| LyndonWordPy { inner })
            .collect()
    }
}

/// Python module providing high-performance Lyndon word operations.
///
/// This module exports the core classes for working with Lyndon words:
/// - `Sort`: Ordering schemes for Lyndon word generation
/// - `LyndonWord`: Individual Lyndon word manipulation
/// - `LyndonBasis`: Factory for generating Lyndon word bases
#[pymodule]
fn lyndon_py(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<SortPy>()?;
    m.add_class::<LyndonBasisPy>()?;
    m.add_class::<LyndonWordPy>()?;
    Ok(())
}
