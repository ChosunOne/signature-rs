use std::{
    collections::HashSet,
    fmt::Display,
    hash::{DefaultHasher, Hash, Hasher},
    ops::{Mul, Neg},
};

use commutator_rs::{Commutator, CommutatorTerm, FormalIndeterminate};
use lyndon_rs::LyndonWord;
use ordered_float::NotNan;
use pyo3::{exceptions::PyValueError, prelude::*, types::PyType};

use crate::formal_indeterminate::FormalIndeterminatePy;

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

impl CommutatorTermPy {
    #[must_use]
    pub fn inner(&self) -> &CommutatorTerm<NotNan<f32>, u8> {
        &self.inner
    }
}

#[pymethods]
impl CommutatorTermPy {
    #[new]
    #[pyo3(signature = (coefficient, atom=None, left=None, right=None))]
    pub fn new(
        coefficient: NotNan<f32>,
        atom: Option<u8>,
        left: Option<Self>,
        right: Option<Self>,
    ) -> PyResult<Self> {
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

    #[classmethod]
    pub fn from_lyndon_word(
        _: Bound<'_, PyType>,
        lyndon_word: &Bound<'_, PyAny>,
    ) -> PyResult<Self> {
        let letters = lyndon_word.getattr("letters")?.extract::<Vec<u8>>()?;
        let lyndon_word = LyndonWord::try_from(letters)
            .map_err(|_| PyValueError::new_err("Invalid lyndon word"))?;
        Ok(Self {
            inner: CommutatorTerm::from(&lyndon_word),
        })
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

    #[must_use]
    pub fn __neg__(&self) -> Self {
        Self {
            inner: self.inner().neg(),
        }
    }

    #[must_use]
    pub fn __hash__(&self) -> u64 {
        let mut hasher = DefaultHasher::new();
        self.inner.hash(&mut hasher);
        hasher.finish()
    }

    #[must_use]
    pub fn degree(&self) -> usize {
        self.inner.degree()
    }

    #[getter]
    #[must_use]
    pub fn get_coefficient(&self) -> f32 {
        **self.inner.coefficient()
    }

    #[setter]
    pub fn set_coefficient(&mut self, coefficient: NotNan<f32>) {
        *self.inner.coefficient_mut() = coefficient;
    }

    #[must_use]
    pub fn commutator(&self, rhs: &Self) -> Self {
        Self {
            inner: self.inner.commutator(&rhs.inner),
        }
    }

    #[getter]
    #[must_use]
    pub fn get_atom(&self) -> Option<u8> {
        self.inner.atom().copied()
    }

    #[setter]
    pub fn set_atom(&mut self, atom: u8) {
        if let Some(x) = self.inner.atom_mut() {
            *x = atom;
        }
    }

    #[getter]
    #[must_use]
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

    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.inner.is_zero()
    }

    #[must_use]
    pub fn unit(&self) -> Self {
        Self {
            inner: self.inner.unit(),
        }
    }

    pub fn lyndon_sort(&mut self) {
        self.inner.lyndon_sort();
    }

    #[must_use]
    pub fn jacobi_identity(&self) -> Option<(Self, Self)> {
        let (a, b) = self.inner.jacobi_identity()?;
        (Self { inner: a }, Self { inner: b }).into()
    }

    #[must_use]
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

    #[must_use]
    pub fn formal_indeterminates(&self) -> Vec<FormalIndeterminatePy> {
        Vec::<FormalIndeterminate<u8, NotNan<f32>>>::from(&self.inner)
            .into_iter()
            .map(|x| FormalIndeterminatePy { inner: x })
            .collect()
    }
}
