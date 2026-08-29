use std::collections::{HashMap, HashSet};
use std::fmt::{Debug, Display};
use std::hash::Hash;
use std::marker::PhantomData;
use std::ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign};
use std::sync::Arc;

use commutator_rs::{Commutator, CommutatorTerm, comm};
use lyndon_rs::generators::Generator;
use lyndon_rs::lyndon::LyndonWord;
use num_traits::{One, Zero};

use crate::feasible_decompositions::{self, FeasibleDecompositions};

#[cfg(feature = "progress")]
use indicatif::{ProgressBar, ProgressStyle};

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

impl<
    T: Clone + Ord + Generator + Hash + Eq,
    U: Clone + One + Zero + Eq + MulAssign + Neg<Output = U> + Hash + AddAssign + Ord,
> LieSeries<T, U>
{
    #[must_use]
    #[cfg_attr(
        feature = "tracing",
        tracing::instrument(skip_all, fields(basis_len = basis.len()))
    )]
    pub fn new(basis: Vec<LyndonWord<T>>, coefficients: Vec<U>) -> Self {
        let mut commutator_basis = Vec::<CommutatorTerm<U, T>>::with_capacity(basis.len());
        let max_degree = if basis.is_empty() {
            0
        } else {
            basis[basis.len() - 1].len()
        };
        for word in &basis {
            commutator_basis.push(CommutatorTerm::from(word));
        }
        #[cfg(feature = "tracing")]
        tracing::debug!(
            max_degree,
            basis_len = commutator_basis.len(),
            "converted Lyndon basis to commutator terms"
        );
        let commutator_basis_set = commutator_basis
            .iter()
            .map(commutator_rs::CommutatorTerm::unit_hash)
            .collect::<HashSet<_>>();

        let mut commutator_basis_index_map = HashMap::new();
        for (i, a) in commutator_basis.iter().enumerate() {
            commutator_basis_index_map.insert(a.unit_hash(), i);
        }
        #[cfg(feature = "tracing")]
        tracing::debug!("built commutator basis maps");

        // Per-word letter contents: the table is grouped by (target degree,
        // content class), which requires each word's letter multiset.
        let mut alphabet: Vec<T> = basis
            .iter()
            .filter(|w| w.len() == 1)
            .map(|w| w.letters[0].clone())
            .collect();
        alphabet.sort();
        alphabet.dedup();
        let basis_contents: Vec<Vec<u8>> = basis
            .iter()
            .map(|w| {
                let mut c = vec![0u8; alphabet.len()];
                for l in &w.letters {
                    let k = alphabet
                        .binary_search(l)
                        .expect("basis letter missing from alphabet");
                    c[k] += 1;
                }
                c
            })
            .collect();
        let mut feasible_builder = feasible_decompositions::Builder::new(&basis_contents);

        #[cfg(feature = "tracing")]
        let _structure_constants_span = tracing::debug_span!(
            "compute_structure_constants",
            pairs = basis.len() * basis.len(),
            max_degree
        )
        .entered();
        #[cfg(feature = "progress")]
        let pb = {
            let style = ProgressStyle::with_template(
                "[{elapsed_precise}] [{bar:35.cyan/blue}] {pos}/{len} structure-constant rows ({msg})",
            )
            .unwrap()
            .progress_chars("=>-");
            let pb = ProgressBar::new(basis.len() as u64).with_style(style);
            pb
        };
        for (i, a) in commutator_basis.iter().enumerate() {
            #[cfg(feature = "progress")]
            {
                pb.set_message(format!("degree {}", a.degree()));
                pb.inc(1);
            }
            #[cfg(feature = "tracing")]
            if i == 0 || commutator_basis[i - 1].degree() != a.degree() {
                tracing::debug!(
                    row = i,
                    degree = a.degree(),
                    "computing structure constants"
                );
            }
            for (j, b) in commutator_basis.iter().enumerate() {
                if j <= i || max_degree < a.degree() + b.degree() {
                    continue;
                }
                #[cfg(feature = "tracing")]
                tracing::trace!(i, j, "computing commutator term");
                let mut term = comm![a, b];
                #[cfg(feature = "tracing")]
                tracing::trace!(i, j, "sorting commutator term");
                term.lyndon_sort();

                // For some non-basis term T, and its decomposition to basis terms A, B, and C, ...
                // T -> [c_1 * A, c_2 * B, c_3 * C, ...]
                #[cfg(feature = "tracing")]
                tracing::trace!(i, j, "decomposing into Lyndon basis");
                let basis_terms = term
                    .lyndon_basis_decomposition(&commutator_basis_set)
                    .into_iter()
                    .filter(|x| !x.coefficient().is_zero())
                    .collect::<Vec<_>>();
                // Get the coefficients [c_1, c_2, c_3, ...]
                let basis_term_coefficients = basis_terms
                    .iter()
                    .map(|x| x.coefficient().clone())
                    .collect::<Vec<_>>();
                // Get the indices i of A, B, C, ...
                // [i_A, i_B, i_C, ...]
                let basis_term_indices = basis_terms
                    .into_iter()
                    .map(|x| commutator_basis_index_map[&x.unit_hash()])
                    .collect::<Vec<_>>();
                #[cfg(feature = "tracing")]
                tracing::trace!(
                    i,
                    j,
                    terms = basis_term_indices.len(),
                    "decomposed commutator into basis terms"
                );
                feasible_builder.push(i, j, &basis_term_indices, &basis_term_coefficients);
            }
        }

        #[cfg(feature = "tracing")]
        tracing::debug!(
            basis_len = basis.len(),
            "finished computing structure constants"
        );
        #[cfg(feature = "progress")]
        pb.finish_with_message("done");
        Self {
            basis: Arc::new(basis),
            commutator_basis: Arc::new(commutator_basis),
            commutator_basis_index_map: Arc::new(commutator_basis_index_map),
            coefficients,
            feasible_decompositions: Arc::new(feasible_builder.finish()),
            max_degree,
        }
    }
}

impl<T, U> LieSeries<T, U> {
    #[must_use]
    pub fn with_coefficients(&self, coefficients: Vec<U>) -> Self {
        Self {
            basis: Arc::clone(&self.basis),
            commutator_basis: Arc::clone(&self.commutator_basis),
            feasible_decompositions: Arc::clone(&self.feasible_decompositions),
            commutator_basis_index_map: Arc::clone(&self.commutator_basis_index_map),
            coefficients,
            max_degree: self.max_degree,
        }
    }

    /// Decomposition of the bracket `[basis[i], basis[j]]` onto the basis,
    /// for pairs stored in the feasible table (`i != j` and degree-feasible).
    /// Returns parallel slices (basis indices, non-zero coefficients).
    #[must_use]
    pub fn decomposition(&self, i: usize, j: usize) -> Option<(&[u32], &[U], bool)> {
        self.feasible_decompositions.get(i, j)
    }
}

/// Minimum visited entries before a batch dispatches to the parallel
/// sweep; below this the fork-join costs more than the work.
const PARALLEL_ENTRY_THRESHOLD: usize = 65536;
/// Target size of a parallel bundle, in visited entries.
const BUNDLE_TARGET_ENTRIES: usize = 2048;

/// The kernel prologue's output: presence bitsets, the active anagram
/// units (shared through a [`GatingCache`]), and the total visited-entry
/// count.
struct KernelGating {
    presence: Vec<u64>,
    words: usize,
    active: Arc<[ActiveSegment]>,
    total_entries: usize,
}

impl KernelGating {
    #[inline]
    fn presences(&self) -> (&[u64], &[u64]) {
        self.presence.split_at(self.words)
    }
}

/// A maximal run of one unit's entries whose smaller endpoint has degree
/// `p` (entries are sorted `(p, q, i, j)` within a unit and `q = target -
/// p` is forced, so runs are contiguous). Gating per run, not per unit:
/// activating a unit for one `p` would otherwise visit every other `p`'s
/// entries, whose presence tests always fail.
#[derive(Clone, Copy)]
struct ActiveSegment {
    /// The run's entry span in the flat entry array, *including* each
    /// entry's flat successor: `entries[span_start .. span_end]` zipped with
    /// itself shifted by one gives entry `k`'s decomposition range
    /// `entries[k].decomp_start .. entries[k + 1].decomp_start`.
    span_start: u32,
    span_end: u32,
    /// The result scatter region `[rs, re)` — the degree-`target` slice.
    /// Decomposition indices are relative to `rs`.
    rs: u32,
    re: u32,
    /// Orientation flags for this run's `p`: `o1` = `(a = i, b = j)` may
    /// contribute, `o2` = the mirrored `(a = j, b = i)` may.
    o1: bool,
    o2: bool,
    /// Whether this is the last active segment of its unit. Segments of one
    /// unit scatter onto the same target words, so the parallel sweep must
    /// keep them in one bundle; it may only cut a bundle after this flag.
    last_of_unit: bool,
}

/// Memo for [`LieSeries::commutator_coefficients_batch_with_cache`]: maps a
/// job's degree-support masks to its active anagram-unit list.
///
/// Unit activity depends only on *which degrees* each side supports — not on
/// coefficient values or the exact support indices — so the prologue reduces
/// reduces each side's non-zero list to a `[u64; 2]` degree mask (covering
/// degrees 0..127, the same range as the table's per-unit `p_mask`) and walks the unit table once per distinct mask pair.
/// In a log-signature fold every DAG node's support lives in a single degree
/// slice (degree homogeneity), so distinct mask pairs are few and nearly
/// every kernel call after the first fold is a memo hit.
///
/// The cache is valid only for the decomposition table of the series whose
/// prologues populated it.
#[derive(Clone, Default)]
pub struct GatingCache {
    /// Open-addressed linear-probe table keyed by the `(a_deg, b_deg)` mask
    /// pair. Distinct pairs per configuration are few (degree-support
    /// signatures of the DAG's nodes), so a fixed-capacity table with a
    /// cheap multiplicative hash beats a full hash map: the lookup runs per
    /// kernel call, and at small grids the SipHash + bucket walk was a
    /// measurable share of the fold.
    slots: Vec<Slot>,
}

/// Key + value for one [`GatingCache`] slot. `None` key = empty; the zero
/// mask pair is a legitimate key, so emptiness is tracked out of band.
#[derive(Clone, Default)]
struct Slot {
    key: Option<([u64; 2], [u64; 2])>,
    value: (Arc<[ActiveSegment]>, usize),
}

impl GatingCache {
    /// Looks up `(a_deg, b_deg)`, returning the memoized
    /// `(active segments, visited entry count)`.
    #[inline]
    fn get(&self, key: ([u64; 2], [u64; 2])) -> Option<&(Arc<[ActiveSegment]>, usize)> {
        let cap = self.slots.len();
        if cap == 0 {
            return None;
        }
        let mut idx = Self::hash(&key) & (cap - 1);
        loop {
            let slot = &self.slots[idx];
            match &slot.key {
                Some(k) if *k == key => return Some(&slot.value),
                Some(_) => idx = (idx + 1) & (cap - 1),
                None => return None,
            }
        }
    }

    /// Inserts, growing the table to the next power of two when full.
    fn insert(&mut self, key: ([u64; 2], [u64; 2]), value: (Arc<[ActiveSegment]>, usize)) {
        if self.slots.len() * 3 / 4 <= self.len() {
            let old = std::mem::take(&mut self.slots);
            self.slots = vec![Slot::default(); (old.len() * 2).max(8)];
            for slot in old {
                if let Some(k) = slot.key {
                    self.insert_unchecked(k, slot.value);
                }
            }
        }
        self.insert_unchecked(key, value);
    }

    fn insert_unchecked(&mut self, key: ([u64; 2], [u64; 2]), value: (Arc<[ActiveSegment]>, usize)) {
        let cap = self.slots.len();
        let mut idx = Self::hash(&key) & (cap - 1);
        while self.slots[idx].key.is_some() {
            idx = (idx + 1) & (cap - 1);
        }
        self.slots[idx] = Slot { key: Some(key), value };
    }

    fn len(&self) -> usize {
        self.slots.iter().filter(|s| s.key.is_some()).count()
    }

    /// Multiplicative fold of the four mask words — cheaper than SipHash
    /// and collision-free enough for a table holding a handful of entries.
    #[inline]
    fn hash(key: &([u64; 2], [u64; 2])) -> usize {
        let mut h = 0xcbf29ce484222325u64;
        for w in [&key.0, &key.1] {
            for x in w {
                h = (h ^ x).wrapping_mul(0x100000001b3);
            }
        }
        h as usize
    }
}

/// Shared handle for the parallel sweep's result writes.
///
/// SAFETY (of the `Send`/`Sync` impls): bundles of anagram units own
/// disjoint index sets — a basis word's content determines the single unit
/// that writes it — so concurrent `scatter_add`s through `ptr` touch
/// disjoint addresses. `U: Send` covers the values moving across threads.
struct RawResult<'a, U> {
    ptr: *mut U,
    _marker: PhantomData<&'a mut [U]>,
}

unsafe impl<U: Send> Send for RawResult<'_, U> {}
unsafe impl<U: Send> Sync for RawResult<'_, U> {}

impl<U> RawResult<'_, U>
where
    U: AddAssign,
{
    #[inline(always)]
    unsafe fn scatter_add(&self, index: usize, value: U) {
        // SAFETY: callers guarantee `index` is in bounds and disjoint across
        // concurrent tasks (the job's result buffer and the unit partition).
        unsafe { *self.ptr.add(index) += value };
    }
}

/// One independent commutation: operand slices plus the destination buffer.
///
/// SAFETY contract for [`LieSeries::commutator_coefficients_batch`]: the
/// `result` buffers of the jobs passed to one batch call must be pairwise
/// disjoint (in a fold these are distinct DAG-node buffers). Within a job,
/// the anagram partition makes the units conflict-free.
pub struct KernelJob<'a, U> {
    /// The left operand's coefficients.
    pub a: &'a [U],
    // SAFETY: `Sync` is sound because the batch's tasks only *read* the
    // operand slices and write through `result` at indices owned by exactly
    // one task (disjoint buffers across jobs, anagram-disjoint units within
    // a job); `U: Send` covers the written values crossing threads.
    //
    // (The unsafe impls follow the struct definition.)
    /// The left operand's non-zero indices (superset of its support).
    pub a_nonzero: &'a [usize],
    /// The right operand's coefficients.
    pub b: &'a [U],
    /// The right operand's non-zero indices.
    pub b_nonzero: &'a [usize],
    /// The destination buffer, as a raw pointer because the batch's parallel
    /// tasks write disjoint regions of it concurrently. Must remain valid
    /// for `result_len` elements for the duration of the call.
    pub result: *mut U,
    /// The destination buffer's length.
    pub result_len: usize,
}

// SAFETY: see the comment inside `KernelJob` — tasks write only indices
// owned by their job, and only values of `U: Send`.
unsafe impl<U: Send> Send for KernelJob<'_, U> {}
unsafe impl<U: Send> Sync for KernelJob<'_, U> {}


/// Batch kernel: evaluates several independent commutation jobs in one
/// parallel dispatch. Jobs must have disjoint `result` buffers; per-word
/// accumulation order matches the serial sweep exactly.
pub fn commutator_coefficients_batch<T, U>(
    a_series: &LieSeries<T, U>,
    jobs: &mut [KernelJob<'_, U>],
) where
    T: Clone + Ord + Generator + Hash,
    U: Clone
        + Default
        + One
        + Zero
        + Eq
        + MulAssign
        + Neg<Output = U>
        + Hash
        + Mul<Output = U>
        + AddAssign
        + Send
        + Sync,
{
    let mut cache = GatingCache::default();
    commutator_coefficients_batch_with_cache(a_series, jobs, &mut cache);
}

/// [`commutator_coefficients_batch`] with a caller-owned [`GatingCache`]
/// that persists across calls, amortizing the gating scan over repeated
/// jobs with equal degree support.
pub fn commutator_coefficients_batch_with_cache<T, U>(
    a_series: &LieSeries<T, U>,
    jobs: &mut [KernelJob<'_, U>],
    cache: &mut GatingCache,
) where
    T: Clone + Ord + Generator + Hash,
    U: Clone
        + Default
        + One
        + Zero
        + Eq
        + MulAssign
        + Neg<Output = U>
        + Hash
        + Mul<Output = U>
        + AddAssign
        + Send
        + Sync,
{
    use rayon::prelude::*;

    // Prologue per job (serial): presence bitsets, degree masks, active
    // units (memoized). Cheap relative to the sweep.
    let gateways: Vec<KernelGating> = jobs
        .iter()
        .map(|j| {
            LieSeries::<T, U>::kernel_prologue_cached(
                a_series,
                j.a_nonzero,
                j.b_nonzero,
                cache,
            )
        })
        .collect();
    let total: usize = gateways.iter().map(|g| g.total_entries).sum();

    if total < PARALLEL_ENTRY_THRESHOLD || rayon::current_num_threads() == 1 {
        for (job, gating) in jobs.iter_mut().zip(&gateways) {
            // SAFETY: the job's result buffer is valid for `result_len`
            // elements (the struct's contract) and is exclusively ours here.
            let result = unsafe { std::slice::from_raw_parts_mut(job.result, job.result_len) };
            LieSeries::sweep_units_serial(a_series, job.a, job.b, gating, result, &mut NoSink);
        }
        return;
    }

    // SAFETY: each job's `result` is a distinct buffer (disjoint across jobs
    // by the caller's contract, disjoint across a job's units by the anagram
    // partition), so the concurrent writes never alias. The buffers are not
    // otherwise accessed during the sweep.
    let writers: Vec<RawResult<U>> = jobs
        .iter()
        .map(|j| RawResult {
            ptr: j.result,
            _marker: PhantomData,
        })
        .collect();

    // Flatten (job, active unit) pairs, carrying each unit's orientation
    // flags, into bundles of roughly `BUNDLE_TARGET_ENTRIES` entries. Units
    // stay in kernel order within a bundle, preserving the per-word
    // accumulation order.
    let table = &a_series.feasible_decompositions;
    // Hoisted flat-table views: every active unit's span indexes directly
    // into these, so neither sweep touches the unit table again.
    let entries_slice = table.entries();
    let decomp_indices_slice = table.decomp_indices_rel();
    let decomp_coeffs_slice = table.decomp_coeffs();
    let mut tasks: Vec<(usize, ActiveSegment)> = Vec::with_capacity(total);
    for (ji, gating) in gateways.iter().enumerate() {
        for au in gating.active.iter() {
            tasks.push((ji, *au));
        }
    }
    let mut bundles: Vec<Vec<(usize, ActiveSegment)>> = vec![Vec::new()];
    let mut cur = 0usize;
    for task in tasks {
        let entries = (task.1.span_end - task.1.span_start - 1) as usize;
        if cur >= BUNDLE_TARGET_ENTRIES && task.1.last_of_unit {
            bundles.push(Vec::new());
            cur = 0;
        }
        bundles.last_mut().unwrap().push(task);
        cur += entries;
    }

    bundles.par_iter().for_each(|bundle| {
        for &(ji, au) in bundle {
            let job = &jobs[ji];
            let writer = &writers[ji];
            let gating = &gateways[ji];
            let (a_present, b_present) = gating.presences();
            let rs = au.rs as usize;
            let span = &entries_slice[au.span_start as usize..au.span_end as usize];
            for (entry, next) in span[..span.len() - 1].iter().zip(span[1..].iter()) {
                let (i, j) = (entry.i as usize, entry.j as usize);
                let p_active = au.o1
                    && a_present[i / 64] & (1u64 << (i % 64)) != 0
                    && b_present[j / 64] & (1u64 << (j % 64)) != 0;
                let q_active = au.o2
                    && a_present[j / 64] & (1u64 << (j % 64)) != 0
                    && b_present[i / 64] & (1u64 << (i % 64)) != 0;
                if !p_active && !q_active {
                    continue;
                }
                let term = if p_active {
                    let mut t = job.a[i].clone() * job.b[j].clone();
                    if q_active {
                        t += -(job.a[j].clone() * job.b[i].clone());
                    }
                    t
                } else {
                    -(job.a[j].clone() * job.b[i].clone())
                };
                let from = entry.decomp_start as usize;
                let to = next.decomp_start as usize;
                for (&rel, c) in decomp_indices_slice[from..to]
                    .iter()
                    .zip(&decomp_coeffs_slice[from..to])
                {
                    // SAFETY: `rs + rel` belongs to this unit's target slice
                    // and to no other unit's (content homogeneity), and is
                    // < job.result_len by the table invariant.
                    unsafe {
                        writer.scatter_add(rs + rel as usize, c.clone() * term.clone());
                    }
                }
            }
        }
    });
}

/// Receiver for the absolute basis indices a kernel call scatters onto.
/// The collecting implementation deduplicates through a dirty bitset.
pub(crate) trait ScatterSink {
    #[inline(always)]
    fn scatter(&mut self, _index: usize) {}
}

pub(crate) struct NoSink;

impl ScatterSink for NoSink {}

pub(crate) struct CollectSink<'a> {
    dirty: &'a mut [u64],
    out: &'a mut Vec<usize>,
    cutoff: usize,
}

impl ScatterSink for CollectSink<'_> {
    #[inline(always)]
    fn scatter(&mut self, index: usize) {
        if index >= self.cutoff {
            return;
        }
        let word = &mut self.dirty[index / 64];
        let bit = 1u64 << (index % 64);
        if *word & bit == 0 {
            *word |= bit;
            self.out.push(index);
        }
    }
}

fn scatter_decomposition<U: Clone + Mul<Output = U> + AddAssign, S: ScatterSink>(
    basis_indices: &[u32],
    basis_coefficients: &[U],
    t: &U,
    result: &mut [U],
    block_start: usize,
    sink: &mut S,
) {
    // SAFETY: `basis_index` is relative to the result slice handed in — for   Not Committed Yet
    // the commutation kernel, the degree-`target` slice of the full result —
    // and comes from the structure-constant table built in `LieSeries::new`,
    // so it is always < result.len() for the duration of the call.
    for (&basis_index, basis_coefficient) in basis_indices.iter().zip(basis_coefficients) {
        sink.scatter(block_start + basis_index as usize);
        *unsafe { result.get_unchecked_mut(basis_index as usize) } +=
            basis_coefficient.clone() * t.clone();
    }
}

impl<
    T: Clone + Ord + Generator + Hash + Eq,
    U: Clone + Default + One + Zero + Eq + MulAssign + Neg<Output = U> + Hash + AddAssign,
> LieSeries<T, U>
{
    /// Indices of non-zero coefficients that the commutation kernel will
    /// actually use: non-zero values on basis elements below `max_degree`.
    pub fn nonzero_coefficient_indices(&self, coefficients: &[U]) -> Vec<usize> {
        let cutoff = self
            .feasible_decompositions
            .degree_start(self.max_degree)
            .min(coefficients.len());
        coefficients[..cutoff]
            .iter()
            .enumerate()
            .filter(|(_, c)| !c.is_zero())
            .map(|x| x.0)
            .collect()
    }

    #[cfg_attr(
        feature = "tracing",
        tracing::instrument(
            level = "debug",
            name = "lie_series_commutator_coefficients",
            skip_all,
            fields(
                basis_len = a_series.basis.len(),
                max_degree = a_series.max_degree
            )
        )
    )]
    pub fn commutator_coefficients(
        a_series: &LieSeries<T, U>,
        a_coefficients: &[U],
        b_coefficients: &[U],
        result_coefficients: &mut [U],
    ) where
        U: Send + Sync,
    {
        let a_nonzero = a_series.nonzero_coefficient_indices(a_coefficients);
        let b_nonzero = a_series.nonzero_coefficient_indices(b_coefficients);

        LieSeries::commutator_coefficients_with_nonzero(
            a_series,
            a_coefficients,
            &a_nonzero,
            b_coefficients,
            &b_nonzero,
            result_coefficients,
        );
    }

    pub fn commutator_coefficients_with_nonzero_collecting(
        a_series: &LieSeries<T, U>,
        a_coefficients: &[U],
        a_nonzero: &[usize],
        b_coefficients: &[U],
        b_nonzero: &[usize],
        result_coefficients: &mut [U],
        dirty: &mut [u64],
        targets: &mut Vec<usize>,
    ) {
        for word in dirty.iter_mut() {
            *word = 0;
        }
        targets.clear();
        let cutoff = a_series
            .feasible_decompositions
            .degree_start(a_series.max_degree);
        let gating = Self::kernel_prologue(a_series, a_nonzero, b_nonzero);
        Self::sweep_units_serial(
            a_series,
            a_coefficients,
            b_coefficients,
            &gating,
            result_coefficients,
            &mut CollectSink {
                dirty,
                out: targets,
                cutoff,
            },
        );
    }

    pub fn commutator_coefficients_with_nonzero(
        a_series: &LieSeries<T, U>,
        a_coefficients: &[U],
        a_nonzero: &[usize],
        b_coefficients: &[U],
        b_nonzero: &[usize],
        result_coefficients: &mut [U],
    ) where
        U: Send + Sync,
    {
        let mut job = KernelJob {
            a: a_coefficients,
            a_nonzero,
            b: b_coefficients,
            b_nonzero,
            result: result_coefficients.as_mut_ptr(),
            result_len: result_coefficients.len(),
        };
        commutator_coefficients_batch(a_series, std::slice::from_mut(&mut job));
        // `job` moved into the slice; its fields are all references/raw
        // pointers — nothing to drop.
    }

    fn kernel_prologue(
        a_series: &LieSeries<T, U>,
        a_nonzero: &[usize],
        b_nonzero: &[usize],
    ) -> KernelGating {
        let mut cache = GatingCache::default();
        Self::kernel_prologue_cached(a_series, a_nonzero, b_nonzero, &mut cache)
    }
    fn kernel_prologue_cached(
        a_series: &LieSeries<T, U>,
        a_nonzero: &[usize],
        b_nonzero: &[usize],
        cache: &mut GatingCache,
    ) -> KernelGating {
        // Bitsets of the non-zero indices of both sides, for O(1) presence
        // checks while walking the units. One allocation backs both:
        // `[0..words]` is A's bitset, `[words..2*words]` is B's.
        let words = a_series.basis.len().div_ceil(64);
        let mut presence = vec![0u64; 2 * words];
        let (a_present, b_present) = presence.split_at_mut(words);
        // One fused pass per side: set the presence bit AND the degree-support
        // mask bit (bit `d` of `a_deg` is set iff A carries a non-zero
        // coefficient on some degree-`d` basis word). Unit gating needs only
        // the masks: a unit with degree pairs `(p, t - p)` contributes
        // through orientation (a = i, b = j) iff some `p` has degree-`p`
        // support in A and degree-`t - p` support in B, and through the
        // mirrored orientation in the transpose case. (Degree-`max_degree`
        // words are filtered from the non-zero lists upstream and never
        // appear as pair endpoints, so only degrees below `max_degree` are
        // ever set.)
        let table = &a_series.feasible_decompositions;
        let mut a_deg = [0u64; 2];
        let mut b_deg = [0u64; 2];
        for &i in a_nonzero {
            a_present[i / 64] |= 1u64 << (i % 64);
            let d = table.degree_of(i);
            debug_assert!(d < 128, "degree masks cover degrees 0..127");
            a_deg[d / 64] |= 1u64 << (d % 64);
        }
        for &j in b_nonzero {
            b_present[j / 64] |= 1u64 << (j % 64);
            let d = table.degree_of(j);
            debug_assert!(d < 128, "degree masks cover degrees 0..127");
            b_deg[d / 64] |= 1u64 << (d % 64);
        }

        // Memoized gating pass: a unit's entries are grouped into contiguous
        // p-runs (entries are sorted `(p, q, i, j)` within the unit and
        // `q = target - p` is forced), and a run is kept only when its own
        // `(p, target - p)` degree pair is supported — unit-level gating
        // would drag every other p's entries through presence tests that
        // always fail. Surviving runs carry their pre-packed spans and
        // orientation flags so neither sweep re-derives anything.
        let cached = cache.get((a_deg, b_deg));
        let (active, total_entries) = match cached {
            Some((active, total)) => (active.clone(), *total),
            None => {
                let entries = table.entries();
                let mut active: Vec<ActiveSegment> = Vec::new();
                let mut total_entries = 0usize;
                for unit in table.units().iter() {
                    let t = unit.target as usize;
                    let (rs, re) = table.degree_range(unit.target);
                    let mut run_start = unit.start;
                    let mut cur_p = u8::MAX;
                    let mut last_active = usize::MAX;
                    // Real entries only: `unit.end` is the trailing
                    // sentinel's slot (its decomp_start closes the last
                    // run's last decomposition range via the +1 span).
                    for ei in unit.start..unit.end {
                        let p = table.degree_of(entries[ei as usize].i as usize) as u8;
                        if p == cur_p {
                            continue;
                        }
                        if cur_p != u8::MAX {
                            total_entries += Self::push_run(
                                &mut active, a_deg, b_deg, cur_p, t, rs, re,
                                run_start, ei, &mut last_active,
                            );
                        }
                        cur_p = p;
                        run_start = ei;
                    }
                    if cur_p != u8::MAX {
                        total_entries += Self::push_run(
                            &mut active, a_deg, b_deg, cur_p, t, rs, re,
                            run_start, unit.end, &mut last_active,
                        );
                    }
                    if last_active != usize::MAX {
                        active[last_active].last_of_unit = true;
                    }
                }
                let value = (Arc::from(active), total_entries);
                cache.insert((a_deg, b_deg), value.clone());
                value
            }
        };

        KernelGating {
            presence,
            words,
            active: active.clone(),
            total_entries,
        }
    }
    fn push_run(
        active: &mut Vec<ActiveSegment>,
        a_deg: [u64; 2],
        b_deg: [u64; 2],
        p: u8,
        t: usize,
        rs: usize,
        re: usize,
        run_start: u32,
        run_end: u32,
        last_active: &mut usize,
    ) -> usize {
        let pu = p as usize;
        let qu = t - pu;
        let o1 = a_deg[pu / 64] >> (pu % 64) & 1 != 0 && b_deg[qu / 64] >> (qu % 64) & 1 != 0;
        let o2 = b_deg[pu / 64] >> (pu % 64) & 1 != 0 && a_deg[qu / 64] >> (qu % 64) & 1 != 0;
        if !o1 && !o2 {
            return 0;
        }
        *last_active = active.len();
        active.push(ActiveSegment {
            span_start: run_start,
            span_end: run_end + 1,
            rs: rs as u32,
            re: re as u32,
            o1,
            o2,
            last_of_unit: false,
        });
        (run_end - run_start) as usize
    }
    fn sweep_units_serial<S: ScatterSink>(
        a_series: &LieSeries<T, U>,
        a_coefficients: &[U],
        b_coefficients: &[U],
        gating: &KernelGating,
        result_coefficients: &mut [U],
        sink: &mut S,
    ) {
        let table = &a_series.feasible_decompositions;
        let (a_present, b_present) = gating.presence.split_at(gating.words);
        let entries = table.entries();
        let decomp_indices = table.decomp_indices_rel();
        let decomp_coefficients = table.decomp_coeffs();
        for au in gating.active.iter() {
            let result = &mut result_coefficients[au.rs as usize..au.re as usize];
            let span = &entries[au.span_start as usize..au.span_end as usize];
            for (entry, next) in span[..span.len() - 1].iter().zip(span[1..].iter()) {
                let (i, j) = (entry.i as usize, entry.j as usize);
                let p_active = au.o1
                    && a_present[i / 64] & (1u64 << (i % 64)) != 0
                    && b_present[j / 64] & (1u64 << (j % 64)) != 0;
                let q_active = au.o2
                    && a_present[j / 64] & (1u64 << (j % 64)) != 0
                    && b_present[i / 64] & (1u64 << (i % 64)) != 0;
                if !p_active && !q_active {
                    continue;
                }
                let term = if p_active {
                    let mut t = a_coefficients[i].clone() * b_coefficients[j].clone();
                    if q_active {
                        t += -(a_coefficients[j].clone() * b_coefficients[i].clone());
                    }
                    t
                } else {
                    // Orientation (a = j, b = i) only: `[basis[j], basis[i]]`
                    // is the negation of the stored decomposition.
                    -(a_coefficients[j].clone() * b_coefficients[i].clone())
                };
                let from = entry.decomp_start as usize;
                let to = next.decomp_start as usize;
                scatter_decomposition(
                    &decomp_indices[from..to],
                    &decomp_coefficients[from..to],
                    &term,
                    result,
                    au.rs as usize,
                    sink,
                );
            }
        }
    }

    pub fn commutator_assign(&mut self, other: &Self)
    where
        U: Send + Sync,
    {
        let mut coefficients = vec![U::default(); self.coefficients.len()];
        LieSeries::commutator_coefficients(
            self,
            &self.coefficients,
            &other.coefficients,
            &mut coefficients,
        );
        self.coefficients = coefficients;
    }
}

impl<
    T: Clone + Ord + Generator + Hash,
    U: Clone + Default + One + Zero + Eq + MulAssign + Neg<Output = U> + Hash + AddAssign + Send + Sync,
> Commutator<&Self> for LieSeries<T, U>
{
    type Output = Self;

    /// Calculates the lie bracket `[A, B]` for a lie series for terms within the commutator basis.
    fn commutator(&self, other: &Self) -> Self::Output {
        let mut coefficients = vec![U::default(); self.coefficients.len()];
        LieSeries::<T, U>::commutator_coefficients(
            self,
            // other,
            &self.coefficients,
            &other.coefficients,
            &mut coefficients,
        );
        self.with_coefficients(coefficients)
    }
}

#[cfg(test)]
mod test {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case(2, 2, vec![1, 2, 3], vec![4, 5, 6], vec![0, 0, -3])]
    #[case(2, 2, vec![3, 2, 1], vec![1, 2, 3], vec![0, 0, 4])]
    #[case(2, 3, vec![1, 2, 3, 4, 5], vec![6, 7, 8, 9, 10], vec![0, 0, -5, -10, 5])]
    #[case(3, 3,
        vec![1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        vec![5, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        vec![0, 0, 0, -7, -14, -7, 0, 0, 0, 0, 0, 0, 0, 0])]
    #[case(3, 3,
        vec![1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        vec![0, 0, 0, -7, -14, -7, 0, 0, 0, 0, 0, 0, 0, 0],
        vec![0, 0, 0, 0, 0, 0, -7, -14, 14, 14, 49, 42, -14, 21])]
    #[case(3, 4,
        vec![
            1, 2, 3, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0
        ],
        vec![
            5, 3, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0
        ],
        vec![
            0, 0, 0, -7, -14, -7, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0
        ],
    )]
    #[case(3, 4,
        vec![
            1, 2, 3, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0
        ],
        vec![
            0, 0, 0, -7, -14, -7, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0
        ],
        vec![
            0, 0, 0, 0, 0, 0, -7, -14,
            14, 14, 49, 42, -14, 21, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
        ],
    )]
    fn test_lie_series_commutation(
        #[case] num_generators: usize,
        #[case] basis_depth: usize,
        #[case] a_coefficients: Vec<i128>,
        #[case] b_coefficients: Vec<i128>,
        #[case] expected_coefficients: Vec<i128>,
    ) {
        use lyndon_rs::lyndon::{LyndonBasis, Sort};
        use num_rational::Ratio;

        let basis = LyndonBasis::<u8>::new(num_generators, Sort::Lexicographical)
            .generate_basis(basis_depth);
        let a_coefficients = a_coefficients
            .into_iter()
            .map(Ratio::<i128>::from_integer)
            .collect::<Vec<_>>();
        let a = LieSeries::new(basis.clone(), a_coefficients);

        let b_coefficients = b_coefficients
            .into_iter()
            .map(Ratio::<i128>::from_integer)
            .collect::<Vec<_>>();
        let b = LieSeries::new(basis.clone(), b_coefficients);
        let expected_coefficients = expected_coefficients
            .into_iter()
            .map(Ratio::<i128>::from_integer)
            .collect::<Vec<_>>();

        let series = comm![a, b];
        assert_eq!(series.coefficients.len(), expected_coefficients.len());
        dbg!(&series.coefficients);
        for (i, c) in series.coefficients.iter().enumerate() {
            assert_eq!(
                *c, expected_coefficients[i],
                "{i}: {c:?} != {:?}",
                expected_coefficients[i]
            );
        }
    }

    #[test]
    fn decompositions_contain_no_zero_coefficients() {
        use lyndon_rs::lyndon::{LyndonBasis, Sort};
        use ordered_float::NotNan;

        for (d, m) in [(2usize, 6usize), (3, 5), (4, 4)] {
            let words = LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
            let series = LieSeries::<u8, NotNan<f64>>::new(words, Vec::new());
            assert!(
                series
                    .feasible_decompositions
                    .iter_entries()
                    .flat_map(|(_, _, _, coefficients)| coefficients)
                    .all(|c| !c.is_zero()),
                "zero coefficient found in decompositions for d={d}, m={m}"
            );
        }
    }

    /// Independently recomputes the decomposition of every feasible canonical
    /// pair and pins the compact table's contents (indices, coefficients,
    /// ordering, feasibility, and canonicality) against it, plus the mirrored
    /// (negated) query path.
    #[test]
    fn feasible_decompositions_match_independent_recomputation() {
        use ordered_float::NotNan;
        use lyndon_rs::lyndon::{LyndonBasis, Sort};

        for (d, m) in [(2usize, 6usize), (3, 5), (4, 4)] {
            let words = LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
            let series = LieSeries::<u8, NotNan<f64>>::new(words, Vec::new());
            let degree = |x: usize| series.commutator_basis[x].degree();
            let index_of = |term: &CommutatorTerm<NotNan<f64>, u8>| {
                series.commutator_basis_index_map[&term.unit_hash()]
            };
            let basis_set: HashSet<_> = series
                .commutator_basis
                .iter()
                .map(CommutatorTerm::unit_hash)
                .collect();

            let mut feasible = 0;
            for unit in series.feasible_decompositions.units() {
                let (p, _q, target) = (unit.p_mask, unit.target as usize, unit.target as usize);
                assert_eq!(p.iter().map(|w| w.count_ones()).sum::<u32>() as usize >= 1, true, "empty unit");
                assert_eq!(
                    unit.gamma.iter().map(|&x| x as usize).sum::<usize>(),
                    target,
                    "unit target degree mismatch"
                );
                let entries: Vec<_> = series
                    .feasible_decompositions
                    .unit_iter(unit)
                    .map(|(e, _, _)| (e.i as usize, e.j as usize))
                    .collect();
                assert!(
                    entries.windows(2).all(|w| w[0] < w[1]),
                    "unit {:?} entries not sorted by (i, j)",
                    unit.gamma
                );

                for (i, j) in entries {
                    assert!(i < j, "non-canonical pair stored");
                    assert!(
                        p[(degree(i) / 64) as usize] >> (degree(i) % 64) & 1 == 1,
                        "pair ({i}, {j}) in a unit of the wrong degrees"
                    );
                    assert_eq!(
                        degree(i) + degree(j),
                        target,
                        "pair ({i}, {j}) in a unit of the wrong target degree"
                    );
                    assert!(degree(i) + degree(j) <= m, "infeasible pair stored");

                    let mut term = comm![
                        &series.commutator_basis[i],
                        &series.commutator_basis[j]
                    ];
                    term.lyndon_sort();
                    let expected: Vec<_> = term
                        .lyndon_basis_decomposition(&basis_set)
                        .into_iter()
                        .filter(|x| !x.coefficient().is_zero())
                        .collect();

                    let expected_indices: Vec<u32> = expected
                        .iter()
                        .map(|x| index_of(x) as u32)
                        .collect();
                    let expected_coefficients: Vec<_> = expected
                        .iter()
                        .map(|x| x.coefficient().clone())
                        .collect();

                    // The canonical query returns the stored data unflagged;
                    // comparing through `decomposition` exercises the
                    // degree-block lookup as well as the contents.
                    let (canonical_indices, canonical_coefficients, swapped) =
                        series.decomposition(i, j).expect("canonical pair query");
                    assert!(!swapped, "(i={i}, j={j}) canonical query flagged");
                    assert_eq!(
                        canonical_indices, &expected_indices[..],
                        "(i={i}, j={j}) indices"
                    );
                    assert_eq!(
                        canonical_coefficients, &expected_coefficients[..],
                        "(i={i}, j={j}) coeffs"
                    );
                    feasible += 1;

                    // The mirrored query returns the same (canonical) data,
                    // flagged so the caller negates it into
                    // [basis[j], basis[i]] orientation.
                    let (mirrored_indices, mirrored_coefficients, swapped) =
                        series.decomposition(j, i).expect("mirrored pair query");
                    assert!(swapped, "(i={i}, j={j}) mirrored query not flagged");
                    assert_eq!(mirrored_indices, &expected_indices[..]);
                    assert_eq!(
                        mirrored_coefficients, &expected_coefficients[..],
                        "(i={i}, j={j}) mirrored data differs from canonical"
                    );
                }
            }
            assert_eq!(series.feasible_decompositions.len(), feasible);
        }
    }

    /// The fused kernel evaluates both bracket orientations of every canonical
    /// pair, so `[A, A]` must vanish exactly — every pair contributes
    /// `c * (a_min * a_max - a_max * a_min) = 0`.
    #[test]
    fn commutator_of_series_with_itself_is_zero() {
        use lyndon_rs::lyndon::{LyndonBasis, Sort};
        use ordered_float::NotNan;

        for (d, m) in [(2usize, 6usize), (3, 5)] {
            let words: Vec<LyndonWord<u8>> =
                LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
            let coefficients: Vec<_> = (0..words.len())
                .map(|i: usize| NotNan::new(((i % 7 + 1) as f64) / 3.0 - 1.1).unwrap())
                .collect();
            let series: LieSeries<u8, NotNan<f64>> = LieSeries::new(words, coefficients);
            let result: LieSeries<u8, NotNan<f64>> = series.commutator(&series);
            assert!(
                result.coefficients.iter().all(|c| c.is_zero()),
                "non-zero coefficient in [A, A] for d={d}, m={m}"
            );
        }
    }

    #[test]
    fn commutator_is_antisymmetric() {
        use lyndon_rs::lyndon::{LyndonBasis, Sort};
        use num_rational::Ratio;

        for (d, m) in [(2usize, 7usize), (3, 6)] {
            let words: Vec<LyndonWord<u8>> =
                LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
            let a_coefficients: Vec<_> = (0..words.len())
                .map(|i: usize| Ratio::from_integer(((i * 13 + 5) % 23) as i128 - 11))
                .collect();
            let b_coefficients: Vec<_> = (0..words.len())
                .map(|i: usize| Ratio::from_integer(((i * 29 + 11) % 19) as i128 - 9))
                .collect();
            let a: LieSeries<u8, Ratio<i128>> = LieSeries::new(words.clone(), a_coefficients);
            let b: LieSeries<u8, Ratio<i128>> = LieSeries::new(words, b_coefficients);
            let ab: LieSeries<u8, Ratio<i128>> = a.commutator(&b);
            let ba: LieSeries<u8, Ratio<i128>> = b.commutator(&a);
            for (x, y) in ab.coefficients.iter().zip(&ba.coefficients) {
                assert_eq!(*x, -*y, "antisymmetry violated for d={d}, m={m}");
            }
        }
    }
}

#[cfg(test)]
mod anagram {
    use super::*;
    use lyndon_rs::lyndon::{LyndonBasis, Sort};
    use ordered_float::NotNan;

    /// The free Lie algebra is multigraded by letter content: the bracket of
    /// two content-homogeneous elements is content-homogeneous of the summed
    /// multiset, and the Lyndon basis refines the grading. So every
    /// decomposition word of `[basis[i], basis[j]]` must carry exactly the
    /// letter multiset of `basis[i]` + `basis[j]` — the invariant that lets
    /// the commutation kernel parallelize over target anagram classes with
    /// disjoint writes.
    #[test]
    fn decompositions_are_content_homogeneous() {
        for (d, m) in [(2usize, 8usize), (3, 8), (4, 8)] {
            let words = LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
            let mut letters: Vec<u8> = words
                .iter()
                .filter(|w| w.len() == 1)
                .map(|w| w.letters[0])
                .collect();
            letters.sort_unstable();
            let content = |w: &LyndonWord<u8>| -> Vec<u8> {
                let mut c = vec![0u8; letters.len()];
                for l in &w.letters {
                    let k = letters.iter().position(|a| a == l).unwrap();
                    c[k] += 1;
                }
                c
            };
            let contents: Vec<Vec<u8>> = words.iter().map(content).collect();

            let series = LieSeries::<u8, NotNan<f64>>::new(words, Vec::new());
            let d_len = series.commutator_basis.len();
            for i in 0..d_len {
                for j in (i + 1)..d_len {
                    if let Some((idx, _, _)) = series.decomposition(i, j) {
                        let mut target = contents[i].clone();
                        for k in 0..letters.len() {
                            target[k] += contents[j][k];
                        }
                        for &x in idx {
                            assert_eq!(
                                contents[x as usize], target,
                                "d={d} m={m}: pair ({i}, {j}) decomposition word {x}                                  violates content homogeneity"
                            );
                        }
                    }
                }
            }
        }
    }
}
