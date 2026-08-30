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

        let mut feasible_builder = Builder::new(basis.len());

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

fn scatter_decomposition<U: Clone + Mul<Output = U> + AddAssign>(
    basis_indices: &[u32],
    basis_coefficients: &[U],
    t: &U,
    result: &mut [U],
) {
    // SAFETY: `basis_index` comes from the structure constant table
    // built in `LieSeries::new` as a decomposition onto this same
    // basis, so it is always < basis.len() == result.len() for the
    // duration of the call.
    for (&basis_index, basis_coefficient) in basis_indices.iter().zip(basis_coefficients) {
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
        coefficients
            .iter()
            .enumerate()
            .filter(|(_, c)| !c.is_zero())
            .filter(|(i, _c)| self.commutator_basis[*i].degree() != self.max_degree)
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

    pub fn commutator_coefficients_with_nonzero(
        a_series: &LieSeries<T, U>,
        a_coefficients: &[U],
        a_nonzero: &[usize],
        b_coefficients: &[U],
        b_nonzero: &[usize],
        result_coefficients: &mut [U],
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

        // Rows are sorted by `j` and the non-zero lists are sorted,
        // so a row of `k` can only match within `j` in `[lo, hi]` of
        // the opposite side's support. `row_window` binary searches that range,
        // keeping the scan proportional to the matching entries.
        let (b_lo, b_hi) = match (b_nonzero.first(), b_nonzero.last()) {
            (Some(&lo), Some(&hi)) => (lo as u32, hi as u32),
            _ => return, // b is identically zero: no pair can contribute
        };
        let (a_lo, a_hi) = match (a_nonzero.first(), a_nonzero.last()) {
            (Some(&lo), Some(&hi)) => (lo as u32, hi as u32),
            _ => return, // a is identically zero: no pair can contribute
        };

        // Sweep 1: canonical rows whose smaller index is in A.
        // The stored decomposition is that of `[basis[k], basis[j]]`
        // with `k < j`, so the fused contribution is `c * (a_k * b_j - a_j * b_k)`
        // orientation (a = k, b = j) is active iff `j` is in B,
        // orientation (a = j, b = k) iff additionally `k` is in B and
        // `j` is in A.
        for &k in a_nonzero {
            if b_hi < k as u32 {
                continue;
            }

            let a_k = &a_coefficients[k];
            let k_in_b = b_present[k / 64] & (1u64 << (k % 64)) != 0;
            let b_k = &b_coefficients[k];

            let (lo, hi) = if k_in_b {
                (b_lo.min(a_lo), b_hi.max(a_hi))
            } else {
                (b_lo, b_hi)
            };

            for (j, basis_indices, basis_coefficients) in
                a_series.feasible_decompositions.row_window(k, lo, hi)
            {
                if b_present[j / 64] & (1u64 << (j % 64)) != 0 {
                    let mut t = a_k.clone() * b_coefficients[j].clone();
                    if k_in_b && a_present[j / 64] & (1u64 << (j % 64)) != 0 {
                        t += -(a_coefficients[j].clone() * b_k.clone());
                    }
                    scatter_decomposition(
                        basis_indices,
                        basis_coefficients,
                        &t,
                        result_coefficients,
                    );
                } else if k_in_b && a_present[j / 64] & (1u64 << (j % 64)) != 0 {
                    let t = -(a_coefficients[j].clone() * b_k.clone());
                    scatter_decomposition(
                        basis_indices,
                        basis_coefficients,
                        &t,
                        result_coefficients,
                    );
                }
            }
        }

        // Sweep 2: canonical rows whose smaller index is in B but not in A.
        // Only orientation (a = j, b = k) can be active there
        for &k in b_nonzero {
            if a_present[k / 64] & (1u64 << (k % 64)) != 0 {
                continue;
            }

            if a_hi <= k as u32 {
                continue;
            }

            let b_k = &b_coefficients[k];
            for (j, basis_indices, basis_coefficients) in
                a_series.feasible_decompositions.row_window(k, a_lo, a_hi)
            {
                if a_present[j / 64] & (1u64 << (j % 64)) != 0 {
                    let t = -(a_coefficients[j].clone() * b_k.clone());
                    scatter_decomposition(
                        basis_indices,
                        basis_coefficients,
                        &t,
                        result_coefficients,
                    );
                }
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
                (0..series.feasible_decompositions.num_rows())
                    .flat_map(|i| series.feasible_decompositions.row(i))
                    .flat_map(|(_, _, coefficients)| coefficients)
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
            let d_len = series.commutator_basis.len();
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
            for i in 0..d_len {
                let row: Vec<_> = series.feasible_decompositions.row(i).collect();
                assert!(row.windows(2).all(|w| w[0].0 < w[1].0), "j not ascending");
                assert!(
                    row.iter().all(|(j, _, _)| *j > i),
                    "non-canonical pair stored"
                );

                for (j, indices, coefficients) in &row {
                    assert!(degree(i) + degree(*j) <= m, "infeasible pair stored");

                    let mut term = comm![&series.commutator_basis[i], &series.commutator_basis[*j]];
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

                    assert_eq!(*indices, &expected_indices[..], "(i={i}, j={j}) indices");
                    assert_eq!(
                        *coefficients,
                        &expected_coefficients[..],
                        "(i={i}, j={j}) coeffs"
                    );
                    feasible += 1;

                    let (mirrored_indices, mirrored_coefficients, swapped) =
                        series.decomposition(*j, i).expect("mirrored pair query");

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
