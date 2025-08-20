use lyndon_rs::{LyndonBasis, LyndonWord, Sort};
use pyo3::prelude::*;

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

#[pyclass(name = "LyndonWord", eq)]
#[derive(PartialEq, Eq)]
struct PyLyndonWord {
    inner: LyndonWord<u8>,
}

#[pymethods]
impl PyLyndonWord {
    #[new]
    fn new(letters: Vec<u8>) -> Self {
        Self {
            inner: LyndonWord { letters },
        }
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
