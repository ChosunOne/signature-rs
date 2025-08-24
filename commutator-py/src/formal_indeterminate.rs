use commutator_rs::FormalIndeterminate;
use ordered_float::NotNan;
use pyo3::{
    exceptions::PyValueError,
    prelude::*,
    types::{PyFloat, PyList, PyType},
};
use std::ops::{Mul, Neg};

use crate::commutator::CommutatorTermPy;

#[pyclass(name = "FormalIndeterminate", eq)]
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct FormalIndeterminatePy {
    pub inner: FormalIndeterminate<u8, NotNan<f32>>,
}

impl FormalIndeterminatePy {
    #[must_use]
    pub fn inner(&self) -> &FormalIndeterminate<u8, NotNan<f32>> {
        &self.inner
    }
}

#[pymethods]
impl FormalIndeterminatePy {
    #[new]
    pub fn new(coefficient: f32, symbols: Vec<u8>) -> PyResult<Self> {
        let coefficient = NotNan::try_from(coefficient)
            .map_err(|_| PyValueError::new_err("Received NaN float value"))?;
        Ok(Self {
            inner: FormalIndeterminate {
                coefficient,
                symbols,
            },
        })
    }

    #[classmethod]
    pub fn from_commutator(
        _: &Bound<'_, PyType>,
        commutator: &Bound<'_, CommutatorTermPy>,
    ) -> PyResult<Vec<Self>> {
        let commutator = commutator.extract::<CommutatorTermPy>()?;

        let indeterminates = Vec::<FormalIndeterminate<u8, NotNan<f32>>>::from(commutator.inner())
            .into_iter()
            .map(|f| Self { inner: f })
            .collect();

        Ok(indeterminates)
    }

    pub fn __mul__<'py>(&self, other: &Bound<'py, PyAny>) -> PyResult<Self> {
        if other.is_instance_of::<PyFloat>() {
            let coefficient = other
                .extract::<f32>()?
                .try_into()
                .map_err(|_| PyValueError::new_err("Received NaN float value"))?;
            return Ok(Self {
                inner: FormalIndeterminate {
                    coefficient,
                    symbols: self.inner.symbols.clone(),
                },
            });
        }

        if other.is_instance_of::<Self>() {
            let rhs = other.extract::<Self>()?;
            return Ok(Self {
                inner: (&self.inner).mul(&rhs.inner),
            });
        }
        Err(PyValueError::new_err(
            "Multiplication is only supported with floats and other formal indeterminates",
        ))
    }

    pub fn __rmul__<'py>(&self, other: &Bound<'py, PyAny>) -> PyResult<Self> {
        if other.is_instance_of::<PyFloat>() {
            let coefficient = other
                .extract::<f32>()?
                .try_into()
                .map_err(|_| PyValueError::new_err("Received NaN float value"))?;
            return Ok(Self {
                inner: FormalIndeterminate {
                    coefficient,
                    symbols: self.inner.symbols.clone(),
                },
            });
        }

        if other.is_instance_of::<Self>() {
            let other = other.extract::<Self>()?;
            return Ok(Self {
                inner: (&other.inner).mul(&self.inner),
            });
        }
        Err(PyValueError::new_err(
            "Multiplication is only supported with floats and other formal indeterminates",
        ))
    }

    #[must_use]
    pub fn __neg__(&self) -> Self {
        Self {
            inner: self.inner().neg(),
        }
    }

    #[must_use]
    #[getter]
    pub fn get_coefficient(&self) -> f32 {
        self.inner.coefficient.into()
    }

    #[setter]
    pub fn set_coefficient(&mut self, coefficient: f32) -> PyResult<()> {
        let coefficient = NotNan::try_from(coefficient)
            .map_err(|_| PyValueError::new_err("Received NaN float value."))?;
        self.inner.coefficient = coefficient;
        Ok(())
    }

    #[getter]
    #[must_use]
    pub fn get_symbols(&self, py: Python<'_>) -> PyResult<Py<PyList>> {
        let result = PyList::new(py, &self.inner.symbols)
            .map(Bound::unbind)?
            .clone_ref(py);
        PyResult::Ok(result)
    }

    #[setter]
    pub fn set_symbols(&mut self, symbols: Bound<'_, PyList>) -> PyResult<()> {
        let symbols = symbols.extract::<Vec<u8>>()?;
        self.inner.symbols = symbols;
        Ok(())
    }
}
