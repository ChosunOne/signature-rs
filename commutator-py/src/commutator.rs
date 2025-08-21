use std::{
    collections::HashSet,
    fmt::Display,
    hash::{DefaultHasher, Hash, Hasher},
    ops::Mul,
};

use commutator_rs::{Commutator, CommutatorTerm};
use ordered_float::NotNan;
use pyo3::{exceptions::PyValueError, prelude::*};

#[pyclass(name = "CommutatorTerm", eq, ord, str)]
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CommutatorTermPy {
    inner: CommutatorTerm<NotNan<f32>, u8>,
}

impl Display for CommutatorTermPy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.inner.fmt(f)
    }
}

#[pymethods]
impl CommutatorTermPy {
    #[new]
    pub fn new(
        coefficient: f32,
        atom: Option<u8>,
        left: Option<Self>,
        right: Option<Self>,
    ) -> PyResult<Self> {
        let coefficient = NotNan::try_from(coefficient)
            .map_err(|_| PyValueError::new_err("Recieved NaN float value."))?;
        if atom.is_some() && (left.is_some() || right.is_some()) {
            return Err(PyValueError::new_err(
                "Either `atom` should be supplied or `left` and `right`",
            ));
        }

        if atom.is_none() && (left.is_none() || right.is_none()) {
            return Err(PyValueError::new_err(
                "Both `left` and `right` should be supplied for `Expression` terms",
            ));
        }

        if let Some(a) = atom {
            return Ok(Self {
                inner: CommutatorTerm::Atom {
                    coefficient,
                    atom: a,
                },
            });
        }

        if let (Some(l), Some(r)) = (left, right) {
            return Ok(Self {
                inner: CommutatorTerm::Expression {
                    coefficient,
                    left: Box::new(l.inner),
                    right: Box::new(r.inner),
                },
            });
        }

        unreachable!()
    }

    pub fn __mul__(&self, other: f32) -> PyResult<Self> {
        let other_float = other
            .try_into()
            .map_err(|_| PyValueError::new_err("Received NaN float value."))?;
        return Ok(Self {
            inner: (&self.inner).mul(other_float),
        });
    }

    pub fn __hash__(&self) -> u64 {
        let mut hasher = DefaultHasher::new();
        self.inner.hash(&mut hasher);
        hasher.finish()
    }

    pub fn degree(&self) -> usize {
        self.inner.degree()
    }

    #[getter]
    pub fn get_coefficient(&self) -> f32 {
        **self.inner.coefficient()
    }

    #[setter]
    pub fn set_coefficient(&mut self, coefficient: f32) -> PyResult<()> {
        *self.inner.coefficient_mut() = NotNan::<f32>::try_from(coefficient)
            .map_err(|_| PyValueError::new_err("Recieved NaN float value"))?;
        Ok(())
    }

    pub fn commutator(&self, rhs: &Self) -> Self {
        Self {
            inner: self.inner.commutator(&rhs.inner),
        }
    }

    #[getter]
    pub fn get_left(&self) -> Option<Self> {
        self.inner.left().map(|l| Self { inner: l.clone() })
    }

    #[setter]
    pub fn set_left(&mut self, mut term: Self) {
        self.inner.left_mut().replace(&mut term.inner);
    }

    #[getter]
    pub fn right(&mut self) -> Option<Self> {
        self.inner.right().map(|r| Self { inner: r.clone() })
    }

    #[setter]
    pub fn set_right(&mut self, mut term: Self) {
        self.inner.right_mut().replace(&mut term.inner);
    }

    pub fn is_zero(&self) -> bool {
        self.inner.is_zero()
    }

    pub fn unit(&self) -> Self {
        Self {
            inner: self.inner.unit(),
        }
    }

    pub fn lyndon_sort(&mut self) {
        self.inner.lyndon_sort();
    }

    pub fn jacobi_identity(&self) -> Option<(Self, Self)> {
        let (a, b) = self.inner.jacobi_identity()?;
        (Self { inner: a }, Self { inner: b }).into()
    }

    pub fn lyndon_basis_decomposition(
        &self,
        lyndon_basis_set: HashSet<CommutatorTermPy>,
    ) -> Vec<Self> {
        let basis_set = lyndon_basis_set.into_iter().map(|x| x.inner).collect();
        self.inner
            .lyndon_basis_decomposition(&basis_set)
            .into_iter()
            .map(|x| Self { inner: x })
            .collect()
    }
}
