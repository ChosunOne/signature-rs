use std::{fmt::Display, hash::Hash, ops::Mul};

use lyndon_rs::{LyndonBasis, LyndonWord, LyndonWordError, Sort};
use pyo3::{exceptions::PyValueError, prelude::*, sync::GILOnceCell, types::PyList};

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
struct PyLyndonWordError(LyndonWordError);

impl From<LyndonWordError> for PyLyndonWordError {
    fn from(value: LyndonWordError) -> Self {
        Self(value)
    }
}

impl From<PyLyndonWordError> for PyErr {
    fn from(value: PyLyndonWordError) -> Self {
        match value.0 {
            LyndonWordError::InvalidWord => PyValueError::new_err("Invalid lyndon word"),
            LyndonWordError::InvalidLetter => PyValueError::new_err("Invalid letter"),
        }
    }
}

#[pyclass(name = "Sort", eq)]
#[derive(Copy, Clone, PartialEq, Eq)]
enum PySort {
    Lexicographical,
    Topological,
}

impl From<PySort> for Sort {
    fn from(value: PySort) -> Self {
        match value {
            PySort::Lexicographical => Self::Lexicographical,
            PySort::Topological => Self::Topological,
        }
    }
}

impl From<Sort> for PySort {
    fn from(value: Sort) -> Self {
        match value {
            Sort::Lexicographical => Self::Lexicographical,
            Sort::Topological => Self::Topological,
        }
    }
}

#[pyclass(name = "LyndonWord", eq, ord, str, frozen, hash)]
#[derive(PartialEq, Eq, Hash, PartialOrd, Ord)]
struct PyLyndonWord {
    inner: LyndonWord<u8>,
}

impl Display for PyLyndonWord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.inner.fmt(f)
    }
}

#[pymethods]
impl PyLyndonWord {
    #[new]
    fn new(letters: Vec<u8>) -> PyResult<Self> {
        Ok(Self {
            inner: LyndonWord::try_from(letters).or(Err(PyValueError::new_err("Invalid word")))?,
        })
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }

    fn __mul__(&self, other: &Self) -> PyResult<Self> {
        Ok(Self {
            inner: (&self.inner)
                .mul(&other.inner)
                .map_err(PyLyndonWordError::from)?,
        })
    }

    #[getter]
    fn get_letters(&self, py: Python<'_>) -> PyResult<Py<PyList>> {
        static LETTERS: GILOnceCell<Py<PyList>> = GILOnceCell::new();
        let result = LETTERS
            .get_or_try_init(py, || {
                PyList::new(py, &self.inner.letters).map(Bound::unbind)
            })?
            .clone_ref(py);
        PyResult::Ok(result)
    }

    fn len(&self) -> usize {
        self.inner.len()
    }

    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    fn goldberg(&self) -> Vec<usize> {
        self.inner.goldberg()
    }

    fn right_factors(&self) -> Vec<Self> {
        self.inner
            .right_factors()
            .into_iter()
            .map(|inner| Self { inner })
            .collect()
    }

    fn factorize(&self) -> (Self, Self) {
        let (a, b) = self.inner.factorize();
        (Self { inner: a }, Self { inner: b })
    }
}

#[pyclass(name = "LyndonBasis")]
struct PyLyndonBasis {
    inner: LyndonBasis<u8>,
}

#[pymethods]
impl PyLyndonBasis {
    #[new]
    fn new(alphabet_size: usize, sort: PySort) -> Self {
        Self {
            inner: LyndonBasis::new(alphabet_size, sort.into()),
        }
    }

    #[getter]
    fn get_alphabet_size(&self) -> usize {
        self.inner.alphabet_size
    }

    #[setter]
    fn set_alphabet_size(&mut self, alphabet_size: usize) {
        self.inner.alphabet_size = alphabet_size;
    }

    #[getter]
    fn get_sort(&self) -> PySort {
        self.inner.sort.into()
    }

    #[setter]
    fn set_sort(&mut self, sort: PySort) {
        self.inner.sort = sort.into();
    }

    fn number_of_words_per_degree(&self, max_degree: usize) -> Vec<usize> {
        self.inner.number_of_words_per_degree(max_degree)
    }

    fn generate_basis(&self, max_length: usize) -> Vec<PyLyndonWord> {
        self.inner
            .generate_basis(max_length)
            .into_iter()
            .map(|inner| PyLyndonWord { inner })
            .collect()
    }
}

/// A Python module implemented in Rust.
#[pymodule]
fn lyndon_py(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySort>()?;
    m.add_class::<PyLyndonBasis>()?;
    m.add_class::<PyLyndonWord>()?;
    Ok(())
}
