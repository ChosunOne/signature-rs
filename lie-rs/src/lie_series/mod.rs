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
pub use self::stages::plan_class_sweep_stages;

mod batch;
pub mod cohort;
mod class_ordered;
mod construction;
mod fold_walk;
mod gating;
mod kernel;
mod raw_ops;
mod stages;

#[cfg(test)]
mod anagram;
#[cfg(test)]
mod obj1_probe;
#[cfg(test)]
mod test;

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

/// The kernel prologue's output: the active units of the feasible-
/// decomposition table — each holding its active-entry ticket range, its
/// exact word set (the result positions it writes), and its contribution
/// (pair) count — plus the total pair count (the planner's work unit).
///
/// The division of parallel work is per unit: the unit word sets
/// PARTITION the sweep's write slots (units are unique by (target degree,
/// content) and every entry producing word `w` lives in `w`'s one unit),
/// so any pack cut between units is race-free, and a word's contributions
/// accumulate inside its one unit in table-entry order — exactly the
/// serial sweep's per-word float summation order, independent of pack
/// cuts, slot counts and claim interleaving.
///
/// Cloning is cheap (all hot fields are `Arc` slices) — the cohort sweep
/// stages clone their planned stage's gateways.
#[derive(Clone)]
struct KernelGating {
    /// Active units in table order, each with its ticket range, word-set
    /// range (into [`Self::unit_words`]) and pair weight. Units whose
    /// active entries all have empty rows are omitted (no work).
    units: Arc<[UnitRange]>,
    /// The flat active-entry ticket list, in table order (a subsequence of
    /// the entry stream — the entry table's natural locality). Each unit's
    /// `ticket_start .. ticket_end` slices this list; a unit's sweep
    /// streams its slice sequentially, computing each term once and adding
    /// its row contributions straight into the result buffer.
    tickets: Arc<[u32]>,
    /// The flat ascending active-word position list: the union of the
    /// units' word sets, sorted (units of one degree are ordered by
    /// content bytes, not position, so the union must be sorted for the
    /// sink/scatter-set walk order to stay byte-identical to the per-word
    /// fan-in's). Exactly the set of positions the sweep writes — the
    /// batch scatter sets are this list.
    unit_words: Arc<[u32]>,
    /// Total active contribution (pair) count: the sweep's work unit (one
    /// multiply-add each), used for pack balancing and the slot policy.
    total_pairs: usize,
}

/// One stage of a class-ordered fold schedule: a kernel sweep level (packs
/// of `(job, unit)` tasks over one level's jobs) or a caller-defined block
/// stage (one atomic claim per block index, dispatched through a shared
/// closure — the gather/accumulate phases of a batched fold). A batch of
/// folds is a linear chain of stages; every slot waits at every stage
/// boundary, so a stage's reads observe exactly the writes of all earlier
/// stages (each completion is published with `Release`, observed with
/// `Acquire`).
pub struct ClassBatchStage<'a, U> {
    /// Sweep stages: `(job index within the stage, unit index in that
    /// job's gating)` tasks in sweep order. Block stages: unused.
    tasks: Arc<[(u32, u32)]>,
    /// Balanced pack ranges into `tasks` (sweep stages) or one block per
    /// pack (block stages).
    packs: Arc<[(usize, usize)]>,
    /// Sweep stages: per-job gating.
    gateways: Vec<KernelGating>,
    /// Sweep stages: the stage's jobs, indexed by the tasks' job ids.
    jobs: &'a [KernelJob<'a, U>],
    /// Block stages: `block(block_index, fold_index)` per claimed task.
    /// Blocks write disjoint ranges (the caller's contract); the stage
    /// counters order all cross-stage dependencies.
    block: Option<&'a (dyn Fn(usize, usize) + Send + Sync)>,
    /// Block stages: which fold of the batch this stage serves (the
    /// block task reads the fold's inputs through it). Sweep stages: 0.
    fold: usize,
    /// Sweep stages: run the 4-lane SoA cohort sweep instead of the scalar
    /// one. A cohort stage's jobs carry SoA-interleaved operand/result
    /// pointers (the four lanes' values of class position `p` live at
    /// `[p*4 .. p*4+4]`, lane `l` at `p*4 + l`), and every other plan
    /// input (tasks, packs, gateways, support lists) is lane-invariant —
    /// see `CommutatorDag::fold_batch_cohort` for how such a plan is
    /// built.
    cohort: bool,
}

/// The per-block task/pack shape shared by every block stage of one
/// batch: `blocks` disjoint-range tasks, one pack per task. Built once
/// per batch, then each block stage clones the `Arc`s (refcount bump —
/// no per-stage allocation).
pub struct BlockShape {
    tasks: Arc<[(u32, u32)]>,
    packs: Arc<[(usize, usize)]>,
}

impl BlockShape {
    /// The block shape over `blocks` disjoint-range tasks.
    pub fn new(blocks: usize) -> Self {
        Self {
            tasks: (0..blocks as u32).map(|i| (i, 0)).collect(),
            packs: (0..blocks).map(|i| (i, i + 1)).collect(),
        }
    }
}

impl<'a, U> ClassBatchStage<'a, U> {
    /// A block stage over `shape`'s disjoint-range tasks: pack `p` runs
    /// `task(p, fold)` exactly once.
    pub fn block_stage(
        shape: &BlockShape,
        task: &'a (dyn Fn(usize, usize) + Send + Sync),
        fold: usize,
    ) -> Self {
        Self {
            tasks: Arc::clone(&shape.tasks),
            packs: Arc::clone(&shape.packs),
            gateways: Vec::new(),
            jobs: &[],
            block: Some(task),
            fold,
            cohort: false,
        }
    }

    /// The cohort (SoA) sweep variant of a planned sweep stage: identical
    /// tasks, packs, gateways and job views (all lane-invariant — see the
    /// `cohort` field), dispatched to the 4-lane SoA kernel instead of the
    /// scalar sweep. The stage's jobs' operand/result pointers must be
    /// SoA-interleaved with lane stride 4 and the plan must have been
    /// built with `cohort = true` — see `plan_class_sweep_stages`.
    pub fn cohort_variant(stage: &Self) -> Self {
        assert!(
            stage.cohort,
            "cohort_variant requires a plan built with cohort = true              (SoA-interleaved buffers)"
        );
        Self {
            tasks: Arc::clone(&stage.tasks),
            packs: Arc::clone(&stage.packs),
            gateways: stage.gateways.clone(),
            jobs: stage.jobs,
            block: None,
            fold: 0,
            cohort: true,
        }
    }
}



