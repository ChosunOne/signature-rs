use std::fmt::Display;
use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};

use commutator_rs::{Commutator, CommutatorTerm};
use lie_rs::LieSeries;
use lyndon_rs::LyndonWord;
use ordered_float::NotNan;
use pyo3::types::PyList;
use pyo3::{exceptions::PyValueError, prelude::*};

use crate::{COMMUTATOR_TERM_CLASS, LYNDON_WORD_CLASS};

#[pyclass(name = "LieSeries", str, sequence)]
#[derive(Clone, Debug)]
pub struct LieSeriesPy {
    pub inner: LieSeries<u8, NotNan<f32>>,
}

impl LieSeriesPy {
    #[must_use]
    pub fn inner(&self) -> &LieSeries<u8, NotNan<f32>> {
        &self.inner
    }
}

impl Display for LieSeriesPy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.inner.fmt(f)
    }
}

#[pymethods]
impl LieSeriesPy {
    #[new]
    pub fn new(basis: Vec<Bound<'_, PyAny>>, coefficients: Vec<NotNan<f32>>) -> PyResult<Self> {
        let basis_words = basis
            .into_iter()
            .map(|w| w.getattr("letters"))
            .collect::<PyResult<Vec<_>>>()?
            .into_iter()
            .map(|l| l.extract::<Vec<u8>>())
            .collect::<PyResult<Vec<_>>>()?
            .into_iter()
            .map(LyndonWord::try_from)
            .collect::<Result<Vec<_>, _>>()
            .map_err(|_| PyValueError::new_err("Invalid lyndon word"))?;
        Ok(Self {
            inner: LieSeries::new(basis_words, coefficients),
        })
    }

    #[must_use]
    pub fn __len__(&self) -> usize {
        self.inner.coefficients.len()
    }

    #[must_use]
    pub fn __getitem__(&self, idx: usize) -> NotNan<f32> {
        self.inner[idx]
    }

    pub fn __setitem__(&mut self, idx: usize, coefficient: NotNan<f32>) {
        self.inner[idx] = coefficient;
    }

    #[must_use]
    pub fn __add__(&self, other: &Self) -> Self {
        Self {
            inner: self.inner().add(other.inner()),
        }
    }

    #[must_use]
    pub fn __radd__(&self, other: &Self) -> Self {
        Self {
            inner: other.inner().add(self.inner()),
        }
    }

    pub fn __iadd__(&mut self, other: &Self) {
        self.inner.add_assign(other.inner());
    }

    #[must_use]
    pub fn __sub__(&self, other: &Self) -> Self {
        Self {
            inner: self.inner().sub(other.inner()),
        }
    }

    #[must_use]
    pub fn __rsub__(&self, other: &Self) -> Self {
        Self {
            inner: other.inner().sub(self.inner()),
        }
    }

    pub fn __isub__(&mut self, other: &Self) {
        self.inner.sub_assign(other.inner());
    }

    #[must_use]
    pub fn __mul__(&self, other: NotNan<f32>) -> Self {
        Self {
            inner: self.inner().mul(other),
        }
    }

    #[must_use]
    pub fn __rmul__(&self, other: NotNan<f32>) -> Self {
        Self {
            inner: self.inner().mul(other),
        }
    }

    pub fn __imul__(&mut self, other: NotNan<f32>) {
        self.inner.mul_assign(other);
    }

    #[getter]
    pub fn get_basis(&self, py: Python<'_>) -> PyResult<Vec<Py<PyAny>>> {
        let lyndon_word_class = LYNDON_WORD_CLASS.get_or_try_init(py, || {
            let module = py.import("lyndon_py")?;
            let class = module.getattr("LyndonWord")?;
            PyResult::Ok(class.unbind())
        })?;
        self.inner
            .basis
            .iter()
            .map(|w| PyList::new(py, &w.letters))
            .collect::<PyResult<Vec<_>>>()?
            .into_iter()
            .map(|l| lyndon_word_class.call1(py, (l,)))
            .collect::<PyResult<Vec<_>>>()
    }

    #[getter]
    pub fn get_commutator_basis(&self, py: Python<'_>) -> PyResult<Vec<Py<PyAny>>> {
        self.inner
            .commutator_basis
            .iter()
            .map(|w| build_commutator_py(py, w))
            .collect::<PyResult<Vec<_>>>()
    }

    #[getter]
    #[must_use]
    pub fn get_coefficients(&self) -> Vec<NotNan<f32>> {
        self.inner.coefficients.clone()
    }

    #[getter]
    #[must_use]
    pub fn max_degree(&self) -> usize {
        self.inner.max_degree
    }

    #[must_use]
    pub fn commutator(&self, other: &Self) -> Self {
        Self {
            inner: self.inner().commutator(other.inner()),
        }
    }
}

fn build_commutator_py(
    py: Python<'_>,
    term: &CommutatorTerm<NotNan<f32>, u8>,
) -> PyResult<Py<PyAny>> {
    let commutator_term_class = COMMUTATOR_TERM_CLASS.get_or_try_init(py, || {
        let module = py.import("commutator_py")?;
        let class = module.getattr("CommutatorTerm")?;
        PyResult::Ok(class.unbind())
    })?;
    match term {
        CommutatorTerm::Atom { coefficient, atom } => {
            commutator_term_class.call1::<(f32, u8)>(py, (coefficient.into_inner(), *atom))
        }
        CommutatorTerm::Expression {
            coefficient,
            left,
            right,
        } => {
            let left = build_commutator_py(py, left)?;
            let right = build_commutator_py(py, right)?;
            commutator_term_class.call1::<(f32, Option<u8>, Py<PyAny>, Py<PyAny>)>(
                py,
                (coefficient.into_inner(), None, left, right),
            )
        }
    }
}
