use std::any::TypeId;
use std::collections::{HashMap, HashSet};
use std::fmt::{Debug, Display};
use std::hash::Hash;
use std::marker::PhantomData;
use std::ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign};
use std::sync::atomic::{AtomicU32, AtomicUsize, Ordering};
use std::sync::{Arc, OnceLock};

use commutator_rs::{Commutator, CommutatorTerm, comm};
use lyndon_rs::generators::Generator;
use lyndon_rs::lyndon::LyndonWord;
use num_traits::{One, Zero};

use crate::feasible_decompositions::{
    self, Entry, FeasibleDecompositions, UnitRange, TICKET_INDEX_MASK, TICKET_P_ACTIVE,
    TICKET_Q_ACTIVE,
};

// Re-exported: the class-contiguous ordering handle behind
// `ClassOrderedCommutation`.
pub use crate::feasible_decompositions::ClassOrder;

pub use self::raw_ops::{raw_add_assign, raw_add_assign_ptr, raw_mul};

pub use self::batch::{
    commutator_coefficients_batch, commutator_coefficients_batch_with_cache,
    commutator_coefficients_class_batch_with_cache, IDENTITY_SHIFTS, KernelJob,
};
pub use self::class_ordered::ClassOrderedCommutation;
pub use self::fold_walk::{
    commutator_coefficients_class_fold_with_cache, planned_sweep_entries, run_class_batch,
    run_class_batch_with_work, work_adaptive_slots,
};
pub use self::gating::GatingCache;
pub use self::stages::{BlockShape, ClassBatchStage, plan_class_sweep_stages};

mod batch;
pub mod cohort;
mod class_ordered;
mod construction;
mod fold_walk;
mod gating;
mod kernel;
mod raw_ops;
mod stages;

/// Represents a formal power series in the free Lie algebra.
///
/// A `LieSeries` is a linear combination of Lyndon words (represented as commutator terms)
/// with coefficients from a ring. This structure provides the foundation for computations
/// involving Baker-Campbell-Hausdorff series and other Lie algebraic operations.
pub struct LieSeries<T, U> {
    /// The Lyndon word basis for the free Lie algebra.
    pub basis: Arc<Vec<LyndonWord<T>>>,
    /// The commutator representation of the Lyndon basis elements.
    pub commutator_basis: Arc<Vec<CommutatorTerm<U, T>>>,
    /// Maps arbitrary commutator terms to their decomposition in the basis.
    /// Invariant: no stored decomposition has a zero coefficient.
    pub(crate) feasible_decompositions: Arc<FeasibleDecompositions<U>>,

    /// Maps basis elements to their index positions for efficient lookup.
    pub commutator_basis_index_map: Arc<HashMap<u64, usize>>,
    /// The coefficients corresponding to each basis element.
    pub coefficients: Vec<U>,
    /// The maximum degree of terms included in this series.
    pub max_degree: usize,
}

impl<T: Debug, U: Debug> Debug for LieSeries<T, U> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("LieSeries")
            .field("basis", &self.basis)
            .field("commutator_basis", &self.commutator_basis)
            .field(
                "feasible_decompositions",
                &self.feasible_decompositions.len(),
            )
            .field(
                "commutator_basis_index_map",
                &self.commutator_basis_index_map,
            )
            .field("coefficients", &self.coefficients)
            .field("max_degree", &self.max_degree)
            .finish()
    }
}

impl<T: Display, U: Display + One + PartialEq> Display for LieSeries<T, U> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for (i, (coefficient, basis_term)) in self
            .coefficients
            .iter()
            .zip(self.commutator_basis.iter())
            .enumerate()
        {
            if i == 0 {
                write!(f, "{coefficient} {basis_term}")?;
                continue;
            }
            write!(f, " + {coefficient} {basis_term}")?;
        }
        Ok(())
    }
}

impl<T: Clone, U: Clone> Clone for LieSeries<T, U> {
    fn clone(&self) -> Self {
        Self {
            basis: self.basis.clone(),
            commutator_basis: Arc::clone(&self.commutator_basis),
            commutator_basis_index_map: Arc::clone(&self.commutator_basis_index_map),
            coefficients: self.coefficients.clone(),
            feasible_decompositions: Arc::clone(&self.feasible_decompositions),
            max_degree: self.max_degree,
        }
    }
}

impl<T, U> Index<usize> for LieSeries<T, U> {
    type Output = U;

    fn index(&self, index: usize) -> &Self::Output {
        &self.coefficients[index]
    }
}

impl<
    T: Clone + Ord + Generator + Hash,
    U: Clone + One + Zero + Eq + MulAssign + Neg<Output = U> + Hash,
> Index<LyndonWord<T>> for LieSeries<T, U>
{
    type Output = U;

    fn index(&self, index: LyndonWord<T>) -> &Self::Output {
        let term = CommutatorTerm::<U, T>::from(&index);
        let i = self.commutator_basis_index_map[&term.unit_hash()];
        &self.coefficients[i]
    }
}

impl<
    T: Clone + Ord + Generator + Hash,
    U: Clone + One + Zero + Eq + MulAssign + Neg<Output = U> + Hash,
> Index<&LyndonWord<T>> for LieSeries<T, U>
{
    type Output = U;

    fn index(&self, index: &LyndonWord<T>) -> &Self::Output {
        let term = CommutatorTerm::<U, T>::from(index);
        let i = self.commutator_basis_index_map[&term.unit_hash()];
        &self.coefficients[i]
    }
}

impl<T, U> IndexMut<usize> for LieSeries<T, U> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.coefficients[index]
    }
}

impl<
    T: Clone + Ord + Generator + Hash,
    U: Clone + One + Zero + Eq + MulAssign + Neg<Output = U> + Hash,
> IndexMut<LyndonWord<T>> for LieSeries<T, U>
{
    fn index_mut(&mut self, index: LyndonWord<T>) -> &mut Self::Output {
        let term = CommutatorTerm::<U, T>::from(&index);
        let i = self.commutator_basis_index_map[&term.unit_hash()];
        &mut self.coefficients[i]
    }
}

impl<
    T: Clone + Ord + Generator + Hash,
    U: Clone + One + Zero + Eq + MulAssign + Neg<Output = U> + Hash,
> IndexMut<&LyndonWord<T>> for LieSeries<T, U>
{
    fn index_mut(&mut self, index: &LyndonWord<T>) -> &mut Self::Output {
        let term = CommutatorTerm::<U, T>::from(index);
        let i = self.commutator_basis_index_map[&term.unit_hash()];
        &mut self.coefficients[i]
    }
}

impl<T: Clone, U: Clone + Add<Output = U>> Add for &LieSeries<T, U> {
    type Output = LieSeries<T, U>;

    fn add(self, rhs: Self) -> Self::Output {
        let mut result = self.clone();
        for i in 0..self.coefficients.len() {
            result.coefficients[i] = self.coefficients[i].clone() + rhs.coefficients[i].clone();
        }
        result
    }
}

impl<T, U: Clone + Add<Output = U>> Add for LieSeries<T, U> {
    type Output = LieSeries<T, U>;

    fn add(mut self, rhs: Self) -> Self::Output {
        for i in 0..self.coefficients.len() {
            self.coefficients[i] = self.coefficients[i].clone() + rhs.coefficients[i].clone();
        }
        self
    }
}

impl<T, U: Clone + AddAssign> AddAssign for LieSeries<T, U> {
    fn add_assign(&mut self, rhs: Self) {
        for i in 0..self.coefficients.len() {
            self.coefficients[i] += rhs.coefficients[i].clone();
        }
    }
}

impl<T, U: Clone + AddAssign> AddAssign<&Self> for LieSeries<T, U> {
    fn add_assign(&mut self, rhs: &Self) {
        for i in 0..self.coefficients.len() {
            self.coefficients[i] += rhs.coefficients[i].clone();
        }
    }
}

impl<T: Clone, U: Clone + Sub<Output = U>> Sub for &LieSeries<T, U> {
    type Output = LieSeries<T, U>;

    fn sub(self, rhs: Self) -> Self::Output {
        let mut result = self.clone();
        for i in 0..self.coefficients.len() {
            result.coefficients[i] = self.coefficients[i].clone() - rhs.coefficients[i].clone();
        }
        result
    }
}

impl<T, U: Clone + Sub<Output = U>> Sub for LieSeries<T, U> {
    type Output = LieSeries<T, U>;

    fn sub(mut self, rhs: Self) -> Self::Output {
        for i in 0..self.coefficients.len() {
            self.coefficients[i] = self.coefficients[i].clone() - rhs.coefficients[i].clone();
        }
        self
    }
}

impl<T, U: Clone + SubAssign> SubAssign for LieSeries<T, U> {
    fn sub_assign(&mut self, rhs: Self) {
        for i in 0..self.coefficients.len() {
            self.coefficients[i] -= rhs.coefficients[i].clone();
        }
    }
}

impl<T, U: Clone + SubAssign> SubAssign<&Self> for LieSeries<T, U> {
    fn sub_assign(&mut self, rhs: &Self) {
        for i in 0..self.coefficients.len() {
            self.coefficients[i] -= rhs.coefficients[i].clone();
        }
    }
}

impl<T, U: Clone + Mul<Output = U> + MulAssign> Mul<U> for LieSeries<T, U> {
    type Output = Self;

    fn mul(mut self, rhs: U) -> Self::Output {
        for i in 0..self.coefficients.len() {
            self.coefficients[i] *= rhs.clone();
        }
        self
    }
}

impl<T: Clone, U: Clone + Mul<Output = U> + MulAssign> Mul<U> for &LieSeries<T, U> {
    type Output = LieSeries<T, U>;

    fn mul(self, rhs: U) -> Self::Output {
        let mut result = self.clone();
        for i in 0..self.coefficients.len() {
            result.coefficients[i] *= rhs.clone();
        }
        result
    }
}

impl<T: Clone, U: Clone + Mul<Output = U> + MulAssign> MulAssign<U> for LieSeries<T, U> {
    fn mul_assign(&mut self, rhs: U) {
        for c in &mut self.coefficients {
            *c *= rhs.clone();
        }
    }
}



