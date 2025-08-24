use std::fmt::Display;
use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};

use commutator_rs::Commutator;
use lie_rs::LieSeries;
use lyndon_rs::LyndonWord;
use ordered_float::{FloatIsNan, NotNan};
use pyo3::{exceptions::PyValueError, prelude::*};

#[pyclass(name = "LieSeries", str, sequence)]
#[derive(Clone, Debug)]
pub struct LieSeriesPy {
    pub inner: LieSeries<u8, NotNan<f32>>,
}

impl LieSeriesPy {
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
    pub fn new(basis: Vec<Bound<'_, PyAny>>, coefficients: Vec<f32>) -> PyResult<Self> {
        let coefficients = coefficients
            .into_iter()
            .map(NotNan::try_from)
            .collect::<Result<Vec<_>, FloatIsNan>>()
            .map_err(|_| PyValueError::new_err("Received NaN float value."))?;
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
    pub fn __getitem__(&self, idx: usize) -> f32 {
        self.inner[idx].into()
    }

    pub fn __setitem__(&mut self, idx: usize, coefficient: f32) -> PyResult<()> {
        let coefficient = NotNan::try_from(coefficient)
            .map_err(|_| PyValueError::new_err("Received NaN float value."))?;
        self.inner[idx] = coefficient;
        PyResult::Ok(())
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

    pub fn __mul__(&self, other: f32) -> PyResult<Self> {
        let other = NotNan::try_from(other)
            .map_err(|_| PyValueError::new_err("Received NaN float value."))?;
        PyResult::Ok(Self {
            inner: self.inner().mul(other),
        })
    }

    pub fn __rmul__(&self, other: f32) -> PyResult<Self> {
        let other = NotNan::try_from(other)
            .map_err(|_| PyValueError::new_err("Received NaN float value."))?;
        PyResult::Ok(Self {
            inner: self.inner().mul(other),
        })
    }

    pub fn __imul__(&mut self, other: f32) -> PyResult<()> {
        let other = NotNan::try_from(other)
            .map_err(|_| PyValueError::new_err("Received NaN float value."))?;
        self.inner.mul_assign(other);
        PyResult::Ok(())
    }

    #[must_use]
    pub fn commutator(&self, other: &Self) -> Self {
        Self {
            inner: self.inner().commutator(other.inner()),
        }
    }
}
