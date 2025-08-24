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

#[pyclass(name = "Sort", eq)]
#[repr(u8)]
#[derive(Copy, Clone, PartialEq, Eq)]
pub enum SortPy {
    Lexicographical = 0,
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

    #[must_use]
    pub fn __mul__(&self, other: &Self) -> PyResult<Self> {
        Ok(Self {
            inner: (&self.inner)
                .mul(&other.inner)
                .map_err(LyndonWordPyError::from)?,
        })
    }

    #[getter]
    #[must_use]
    pub fn get_letters(&self, py: Python<'_>) -> PyResult<Py<PyList>> {
        let result = PyList::new(py, &self.inner.letters)
            .map(Bound::unbind)?
            .clone_ref(py);
        PyResult::Ok(result)
    }

    #[must_use]
    pub fn __len__(&self) -> usize {
        self.inner.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    #[must_use]
    pub fn goldberg(&self) -> Vec<usize> {
        self.inner.goldberg()
    }

    #[must_use]
    pub fn right_factors(&self) -> Vec<Self> {
        self.inner
            .right_factors()
            .into_iter()
            .map(|inner| Self { inner })
            .collect()
    }

    #[must_use]
    pub fn factorize(&self) -> (Self, Self) {
        let (a, b) = self.inner.factorize();
        (Self { inner: a }, Self { inner: b })
    }
}

#[pyclass(name = "LyndonBasis")]
#[derive(Clone, Debug, Copy, Default)]
pub struct LyndonBasisPy {
    inner: LyndonBasis<u8>,
}

#[pymethods]
impl LyndonBasisPy {
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

    #[must_use]
    pub fn number_of_words_per_degree(&self, max_degree: usize) -> Vec<usize> {
        self.inner.number_of_words_per_degree(max_degree)
    }

    #[must_use]
    pub fn generate_basis(&self, max_length: usize) -> Vec<LyndonWordPy> {
        self.inner
            .generate_basis(max_length)
            .into_iter()
            .map(|inner| LyndonWordPy { inner })
            .collect()
    }
}

/// A Python module implemented in Rust.
#[pymodule]
fn lyndon_py(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<SortPy>()?;
    m.add_class::<LyndonBasisPy>()?;
    m.add_class::<LyndonWordPy>()?;
    Ok(())
}
