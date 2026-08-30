use std::collections::{HashMap, HashSet};
use std::fmt::{Debug, Display};
use std::hash::Hash;
use std::ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign};
use std::sync::Arc;

use commutator_rs::{Commutator, CommutatorTerm, comm};
use lyndon_rs::generators::Generator;
use lyndon_rs::lyndon::LyndonWord;
use num_traits::{One, Zero};

use crate::feasible_decompositions::{Builder, FeasibleDecompositions};

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

        let basis_degrees = commutator_basis
            .iter()
            .map(|w| u8::try_from(w.degree()).expect("word degree exceeds u8"))
            .collect::<Vec<_>>();
        let mut feasible_builder = Builder::new(&basis_degrees);

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
    ) {
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
        Self::commutation_kernel(
            a_series,
            a_coefficients,
            a_nonzero,
            b_coefficients,
            b_nonzero,
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
    ) {
        Self::commutation_kernel(
            a_series,
            a_coefficients,
            a_nonzero,
            b_coefficients,
            b_nonzero,
            result_coefficients,
            &mut NoSink,
        );
    }

    fn commutation_kernel<S: ScatterSink>(
        a_series: &LieSeries<T, U>,
        a_coefficients: &[U],
        a_nonzero: &[usize],
        b_coefficients: &[U],
        b_nonzero: &[usize],
        result_coefficients: &mut [U],
        sink: &mut S,
    ) {
        #[cfg(feature = "tracing")]
        tracing::debug!(
            nonzero_a = a_nonzero.len(),
            nonzero_b = b_nonzero.len(),
            "computing commutator coefficients"
        );

        let words = a_series.basis.len().div_ceil(64);
        let mut presence = vec![0u64; 2 * words];
        let (a_present, b_present) = presence.split_at_mut(words);

        for &i in a_nonzero {
            a_present[i / 64] |= 1u64 << (i % 64);
        }
        for &j in b_nonzero {
            b_present[j / 64] |= 1u64 << (j % 64);
        }

        let table = &a_series.feasible_decompositions;
        let max_degree = table.max_degree();
        let mut a_deg = vec![false; max_degree + 1];
        let mut b_deg = vec![false; max_degree + 1];
        for &i in a_nonzero {
            a_deg[table.degree_of(i)] = true;
        }
        for &j in b_nonzero {
            b_deg[table.degree_of(j)] = true;
        }

        for block in table.blocks() {
            let (p, q) = (block.p as usize, block.q as usize);
            let o1 = a_deg[p] && b_deg[q];
            let o2 = a_deg[q] && b_deg[p];
            if !o1 && !o2 {
                continue;
            }
            let (rs, re) = table.degree_range(block.target);
            let result = &mut result_coefficients[rs..re];
            let span = table.entry_span(block);
            let decomp_indices = table.decomp_indices_rel();
            let decomp_coefficients = table.decomp_coeffs();
            for (entry, next) in span[..span.len() - 1].iter().zip(span[1..].iter()) {
                let (i, j) = (entry.i as usize, entry.j as usize);
                let p_active = o1
                    && a_present[i / 64] & (1u64 << (i % 64)) != 0
                    && b_present[j / 64] & (1u64 << (j % 64)) != 0;
                let q_active = o2
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
                    -(a_coefficients[j].clone() * b_coefficients[i].clone())
                };
                let from = entry.decomp_start as usize;
                let to = next.decomp_start as usize;
                scatter_decomposition(
                    &decomp_indices[from..to],
                    &decomp_coefficients[from..to],
                    &term,
                    result,
                    rs,
                    sink,
                );
            }
        }
    }

    pub fn commutator_assign(&mut self, other: &Self) {
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
    U: Clone + Default + One + Zero + Eq + MulAssign + Neg<Output = U> + Hash + AddAssign,
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

    #[test]
    fn feasible_decompositions_match_independent_recomputation() {
        use lyndon_rs::lyndon::{LyndonBasis, Sort};
        use ordered_float::NotNan;

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
            for block in series.feasible_decompositions.blocks() {
                let (p, q, target) = (block.p as usize, block.q as usize, block.target as usize);
                assert!(p <= q, "block degrees are not canonical: ({p}, {q})");
                assert_eq!(p + q, target, "block target degree mismatch");
                let entries = series
                    .feasible_decompositions
                    .block_iter(block)
                    .map(|(e, _, _)| (e.i as usize, e.j as usize))
                    .collect::<Vec<_>>();
                assert!(
                    entries.windows(2).all(|w| w[0] < w[1]),
                    "block ({p}, {q}) entries are not sorted by (i, j)"
                );

                for (i, j) in entries {
                    assert!(i < j, "non-canonical pair stored");
                    assert_eq!(
                        series.feasible_decompositions.degree_of(i),
                        p,
                        "pair ({i}, {j}) in a block of the wrong degree"
                    );
                    assert_eq!(
                        series.feasible_decompositions.degree_of(j),
                        q,
                        "pair ({i}, {j}) in a block of the wrong degree"
                    );
                    assert!(degree(i) + degree(j) <= m, "infeasible pair stored");

                    let mut term = comm![&series.commutator_basis[i], &series.commutator_basis[j]];
                    term.lyndon_sort();
                    let expected: Vec<_> = term
                        .lyndon_basis_decomposition(&basis_set)
                        .into_iter()
                        .filter(|x| !x.coefficient().is_zero())
                        .collect();
                    let expected_indices: Vec<u32> =
                        expected.iter().map(|x| index_of(x) as u32).collect();
                    let expected_coefficients: Vec<_> =
                        expected.iter().map(|x| x.coefficient().clone()).collect();

                    let (canonical_indices, canonical_coefficients, swapped) =
                        series.decomposition(i, j).expect("canonical pair query");
                    assert!(!swapped, "(i={i}, j={j}) canonical query flagged");
                    assert_eq!(
                        canonical_indices,
                        &expected_indices[..],
                        "(i={i}, j={j}) indices"
                    );
                    assert_eq!(
                        canonical_coefficients,
                        &expected_coefficients[..],
                        "(i={i}, j={j}) coeffs"
                    );
                    feasible += 1;

                    let (mirrored_indices, mirrored_coefficients, swapped) =
                        series.decomposition(j, i).expect("mirrored pair query");

                    assert!(swapped, "(i={i}, j={j}) mirrored query not flagged");
                    assert_eq!(mirrored_indices, &expected_indices[..]);
                    assert_eq!(
                        mirrored_coefficients,
                        &expected_coefficients[..],
                        "(i={i}, j={j}) mirrored data differs from canonical"
                    );
                }
            }
            assert_eq!(series.feasible_decompositions.len(), feasible);
        }
    }

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
