use lyndon_rs::{
    generators::Generator,
    lyndon::{LyndonBasis, LyndonWord},
};
use num_traits::{FromPrimitive, One, Zero};
use rayon::prelude::*;
use std::{
    collections::HashMap,
    fmt::{Debug, Display},
    hash::Hash,
    marker::PhantomData,
    ops::{AddAssign, Div, Mul, MulAssign, Neg, SubAssign},
};

#[cfg(feature = "progress")]
mod progress_imports {
    pub use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
}

#[cfg(feature = "progress")]
pub use progress_imports::*;

use crate::{
    LieSeries, LieSeriesGenerator,
    bch::{bch_denominator, goldberg_coeff_numerator},
    binomial,
    constants::FACTORIALS,
};

/// Generator for Baker-Campbell-Hausdorff (BCH) series using Lyndon words.
///
/// This structure efficiently computes BCH series coefficients by organizing
/// the Lyndon basis and precomputing factorizations. The BCH formula gives
/// a way to express log(exp(X)exp(Y)) as a series of nested commutators.
pub struct BchSeriesGenerator<T> {
    /// The size of the generator alphabet.
    pub alphabet_size: usize,
    /// Maximum degree of terms to include in the series.
    pub max_degree: usize,
    /// The Lyndon words forming the basis of the free Lie algebra.
    pub basis: Vec<LyndonWord<T>>,
    /// Index of the left factor for each word in the factorization.
    pub left_factor: Vec<usize>,
    /// Index of the right factor for each word in the factorization.
    pub right_factor: Vec<usize>,
    /// Length (degree) of each word in the basis.
    pub word_lengths: Vec<usize>,
    /// Starting indices for words of each degree in the basis.
    pub index_of_degree: Vec<usize>,
    /// Multi-degree information for each basis word.
    pub multi_degree: Vec<usize>,
}

impl<T: Debug> Debug for BchSeriesGenerator<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BchSeriesGenerator")
            .field("alphabet_size", &self.alphabet_size)
            .field("basis", &self.basis)
            .field("max_degree", &self.max_degree)
            .field("basis", &self.basis)
            .field("left_factor", &self.left_factor)
            .field("right_factor", &self.right_factor)
            .field("word_lengths", &self.word_lengths)
            .field("index_of_degree", &self.index_of_degree)
            .field("multi_degree", &self.multi_degree)
            .finish()
    }
}

impl<T: Clone> Clone for BchSeriesGenerator<T> {
    fn clone(&self) -> Self {
        Self {
            alphabet_size: self.alphabet_size,
            max_degree: self.max_degree,
            basis: self.basis.clone(),
            left_factor: self.left_factor.clone(),
            right_factor: self.right_factor.clone(),
            word_lengths: self.word_lengths.clone(),
            index_of_degree: self.index_of_degree.clone(),
            multi_degree: self.multi_degree.clone(),
        }
    }
}

impl<T: Clone + Eq + Hash + Ord + Generator + Send + Sync> BchSeriesGenerator<T> {
    /// Creates a new BCH series generator from a Lyndon basis.
    ///
    /// This constructor precomputes all necessary factorizations and indices
    /// for efficient BCH coefficient computation.
    #[must_use]
    pub fn new(basis: LyndonBasis<T>, max_degree: usize) -> Self {
        #[cfg(feature = "progress")]
        let style = ProgressStyle::with_template(
            "[{eta_precise}] [{bar:35.green/white}] {pos:>2}/{len:2} {msg}",
        )
        .unwrap()
        .progress_chars("=>-");

        let alphabet = T::alphabet(basis.alphabet_size);
        let alphabet_size = basis.alphabet_size;
        let number_of_words_per_degree = basis.number_of_words_per_degree(max_degree);
        let basis = basis.generate_basis(max_degree);

        let mut left_factor = vec![0; basis.len()];
        for (i, lf) in left_factor.iter_mut().enumerate().take(alphabet_size) {
            *lf = i;
        }
        let mut right_factor = vec![0; basis.len()];
        let mut word_lengths = vec![0; basis.len()];
        let mut index_of_degree = vec![0; max_degree + 1];

        for n in 1..=max_degree {
            index_of_degree[n] = index_of_degree[n - 1] + number_of_words_per_degree[n - 1];
        }

        #[cfg(feature = "progress")]
        let pb = ProgressBar::new(basis.len() as u64).with_style(style.clone());
        #[cfg(feature = "progress")]
        {
            pb.set_message("Factorizing Lyndon words ");
            pb.tick();
        }

        basis
            .par_iter()
            .zip(left_factor.par_iter_mut())
            .zip(right_factor.par_iter_mut())
            .zip(word_lengths.par_iter_mut())
            .for_each(|(((word, lf), rf), wl)| {
                if word.len() > 1 {
                    let (l, r) = word.factorize();
                    for (k, basis_k) in basis
                        .iter()
                        .enumerate()
                        .take((index_of_degree[l.len()] - 1) + 1)
                        .skip(index_of_degree[l.len() - 1])
                    {
                        if *basis_k == l {
                            *lf = k;
                            break;
                        }
                    }
                    for (k, basis_k) in basis
                        .iter()
                        .enumerate()
                        .take((index_of_degree[r.len()] - 1) + 1)
                        .skip(index_of_degree[r.len() - 1])
                    {
                        if *basis_k == r {
                            *rf = k;
                            break;
                        }
                    }
                }
                *wl = word.len();
            });

        #[cfg(feature = "progress")]
        pb.finish();

        let mut multi_degree = vec![0; basis.len()];
        let mut alphabet_index = HashMap::new();
        for (i, letter) in alphabet.iter().enumerate() {
            alphabet_index.insert(letter, i);
        }

        #[cfg(feature = "progress")]
        {
            pb.set_length(0);
            pb.set_message("Calculating multi-degree groups ");
            pb.tick();
        }

        for i in 0..basis.len() {
            let mut letter_counts = vec![0usize; alphabet_size];
            let word = &basis[i];
            for j in 0..word_lengths[i] {
                letter_counts[alphabet_index[&word.letters[j]]] += 1;
            }
            multi_degree[i] = tuple_index(&letter_counts);
            #[cfg(feature = "progress")]
            pb.inc(1);
        }
        #[cfg(feature = "progress")]
        pb.finish();

        Self {
            alphabet_size,
            max_degree,
            basis,
            left_factor,
            right_factor,
            word_lengths,
            index_of_degree,
            multi_degree,
        }
    }

    #[must_use]
    pub fn generate_goldberg_coefficient_numerators<
        U: Clone
            + FromPrimitive
            + One
            + Zero
            + Div<Output = U>
            + AddAssign
            + MulAssign
            + Neg<Output = U>
            + PartialEq
            + Send
            + Sync,
    >(
        &self,
    ) -> Vec<U> {
        #[cfg(feature = "progress")]
        let style = ProgressStyle::with_template(
            "[{eta_precise}] [{bar:35.green/white}] {pos:>2}/{len:2} {msg}",
        )
        .unwrap()
        .progress_chars("=>-");

        let mut goldberg_coefficient_numerators = vec![U::zero(); self.basis.len()];

        let alphabet = T::alphabet(1);
        #[cfg(feature = "progress")]
        let pb = ProgressBar::new(self.basis.len() as u64).with_style(style.clone());
        #[cfg(feature = "progress")]
        {
            pb.set_message("Calculating Goldberg Coefficients ");
            pb.tick();
        }

        self.basis
            .par_iter()
            .zip(goldberg_coefficient_numerators.par_iter_mut())
            .for_each(|(word, numer)| {
                let goldberg_partition = word.goldberg();
                let a_first = word.letters[0] == alphabet[0];
                *numer = goldberg_coeff_numerator(&goldberg_partition, a_first);
            });

        #[cfg(feature = "progress")]
        pb.finish();

        goldberg_coefficient_numerators
    }

    fn right_factors(&self, i: usize, kmax: usize) -> (usize, Vec<usize>) {
        let mut j = vec![0; kmax];
        let mut k = 0;
        j[0] = i;
        let mut l = i;
        while k < kmax && self.left_factor[l] == 0 {
            k += 1;
            l = self.right_factor[l];
            j[k] = l;
        }
        (k, j)
    }

    #[allow(clippy::too_many_lines)]
    fn generate_bch_coefficients<
        U: Clone
            + AddAssign
            + Mul<Output = U>
            + MulAssign
            + Neg<Output = U>
            + Div<Output = U>
            + FromPrimitive
            + PartialEq
            + One
            + Zero
            + SubAssign
            + Send
            + Sync,
    >(
        &self,
    ) -> Vec<U> {
        #[cfg(feature = "progress")]
        let style = ProgressStyle::with_template(
            "[{elapsed_precise}] [{bar:35.green/white}] {pos:>2}/{len:2} {msg}",
        )
        .unwrap()
        .progress_chars("=>-");
        #[cfg(feature = "progress")]
        let style_mb = ProgressStyle::with_template(
            "[{eta_precise}] [{bar:35.green/white}] {pos:>2}/{len:2} {msg}",
        )
        .unwrap()
        .progress_chars("=>-");
        let mut bch_coefficient_numerators = self.generate_goldberg_coefficient_numerators::<U>();
        let alphabet = T::alphabet(1);

        // Loop over the words of max degree
        let start = self.index_of_degree[self.max_degree - 1];
        let end = self.index_of_degree[self.max_degree] - 1;
        let h1 = self.multi_degree[start];
        let h2 = self.multi_degree[end];

        #[cfg(feature = "progress")]
        let mb = MultiProgress::new();
        #[cfg(feature = "progress")]
        mb.println("Calculating BCH Coefficients").unwrap();
        #[cfg(feature = "progress")]
        let p1 = mb.add(ProgressBar::new((h2 - h1 + 1) as u64).with_style(style.clone()));
        #[cfg(feature = "progress")]
        let p2 = mb.add(ProgressBar::new(0).with_style(style_mb.clone()));

        #[cfg(feature = "progress")]
        p1.set_message("Calculating by multi-degree ");

        let coeffs_addr = bch_coefficient_numerators.as_mut_ptr() as usize;
        let len = bch_coefficient_numerators.len();

        (h1..=h2).into_par_iter().for_each(|h| {
            let words_in_group = self.basis[start..=end]
                .iter()
                .enumerate()
                .filter(|&(i, _)| self.multi_degree[i + start] == h)
                .map(|(i, _)| i + start)
                .collect::<Vec<_>>();

            let mut matrix_tree = MatrixTree::<T>::new(self.max_degree as u8, self.alphabet_size);
            let mut group_matrix_indices = Vec::with_capacity(words_in_group.len());

            for &j in &words_in_group {
                let (num_factors, right_factors) = self.right_factors(j, self.max_degree);

                let rightmost_factor = right_factors[num_factors];

                group_matrix_indices.push(matrix_tree.append(
                    rightmost_factor,
                    num_factors as u8,
                    &self.left_factor,
                    &self.right_factor,
                    &self.word_lengths,
                ));
            }

            let mut x = vec![0_i64; matrix_tree.matrices.len()];
            let mut stop = 0;
            for idx in 0..words_in_group.len() {
                let i = words_in_group[idx];
                stop = if group_matrix_indices[idx] > stop {
                    group_matrix_indices[idx]
                } else {
                    stop
                };

                let (num_right_factors, right_factors) = self.right_factors(i, self.max_degree);
                let word = &self.basis[i];
                matrix_tree.run(&mut x, word, stop);

                let mut index_of_non_initial_generator = 0;
                while word.letters[index_of_non_initial_generator] == alphabet[0] {
                    index_of_non_initial_generator += 1;
                }

                for y in 0..idx {
                    let j = words_in_group[y];
                    let (previous_num_right_factors, previous_right_factors) =
                        self.right_factors(j, self.max_degree);
                    if index_of_non_initial_generator >= previous_num_right_factors {
                        let d = x[group_matrix_indices[y]];
                        if d != 0 {
                            for k in 0..=previous_num_right_factors.min(num_right_factors) {
                                let previous_word_coefficient =
                                    bch_coefficient_numerators[previous_right_factors[k]].clone();
                                // SAFETY: Different mutli-degree groups
                                // are disjoint, so this is safe.
                                unsafe {
                                    let ptr = coeffs_addr as *mut U;
                                    let coeffs = std::slice::from_raw_parts_mut(ptr, len);
                                    coeffs[right_factors[k]] -= U::from_i64(d)
                                        .expect("Failed to convert from i64")
                                        * previous_word_coefficient;
                                }
                            }
                        }
                    }
                }
            }
        });

        #[cfg(feature = "progress")]
        {
            p1.finish();
            p2.finish();
        }

        bch_coefficient_numerators
            .into_iter()
            .zip(self.basis.iter())
            .map(|(numerator, word)| {
                let n = word.goldberg().iter().sum();
                let denominator = U::from_i128(FACTORIALS[n]).expect("Failed to convert from i128")
                    * bch_denominator::<U>(n);
                numerator / denominator
            })
            .collect()
    }
}

impl<
    T: Clone + Eq + Hash + Ord + Generator + Send + Sync,
    U: Clone
        + AddAssign
        + Div<Output = U>
        + FromPrimitive
        + Neg<Output = U>
        + One
        + Ord
        + Hash
        + Send
        + Sync
        + Zero
        + SubAssign
        + MulAssign,
> LieSeriesGenerator<T, U> for BchSeriesGenerator<T>
{
    fn generate_lie_series(&self) -> LieSeries<T, U> {
        let basis = self.basis.clone();
        let coefficients = self.generate_bch_coefficients::<U>();

        LieSeries::new(basis, coefficients)
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub struct Matrix2x2 {
    a11: u32,
    a12: u32,
    a21: u32,
    a22: u32,
}

pub struct MatrixTree<T: Generator> {
    alphabet_size: usize,
    /// A vector of matrices
    matrices: Vec<Matrix2x2>,
    /// The maximum degree of the words
    degree: u8,
    /// Maps keys to indices in `matrices`.
    memoization_cache: HashMap<usize, usize>,
    _phantom: PhantomData<T>,
}

impl<T: Generator + Eq> MatrixTree<T> {
    #[must_use]
    pub fn new(degree: u8, alphabet_size: usize) -> Self {
        let mut matrices = Vec::with_capacity(alphabet_size * degree as usize);
        let mut memoization_cache = HashMap::new();

        for i in 0..alphabet_size {
            for j in 0..degree {
                let matrix = Matrix2x2::default();
                let key = (i << 8) | j as usize;
                memoization_cache.insert(key, matrices.len());
                matrices.push(matrix);
            }
        }

        Self {
            alphabet_size,
            matrices,
            degree,
            memoization_cache,
            _phantom: PhantomData,
        }
    }

    #[allow(clippy::cast_possible_truncation)]
    pub fn append(
        &mut self,
        word_index: usize,
        l: u8,
        left_factor: &[usize],
        right_factor: &[usize],
        word_lengths: &[usize],
    ) -> usize {
        let key = (word_index << 8) | l as usize;
        if let Some(&v) = self.memoization_cache.get(&key) {
            return v;
        }

        // Recursively decompose the word into its factors
        let a11 = self.append(
            left_factor[word_index],
            l,
            left_factor,
            right_factor,
            word_lengths,
        ) as u32;
        let a12 = self.append(
            left_factor[word_index],
            l + word_lengths[right_factor[word_index]] as u8,
            left_factor,
            right_factor,
            word_lengths,
        ) as u32;
        let a21 = self.append(
            right_factor[word_index],
            l,
            left_factor,
            right_factor,
            word_lengths,
        ) as u32;
        let a22 = self.append(
            right_factor[word_index],
            l + word_lengths[left_factor[word_index]] as u8,
            left_factor,
            right_factor,
            word_lengths,
        ) as u32;

        self.memoization_cache.insert(key, self.matrices.len());
        self.matrices.push(Matrix2x2 { a11, a12, a21, a22 });
        self.matrices.len() - 1
    }

    pub fn run(&self, x: &mut [i64], word: &LyndonWord<T>, mut stop: usize) {
        let alphabet = T::alphabet(self.alphabet_size);
        if stop >= self.matrices.len() {
            stop = self.matrices.len() - 1;
        }
        for k in 0..self.alphabet_size {
            for i in 0..self.degree {
                x[k * self.degree as usize + i as usize] =
                    i64::from(word.letters[i as usize] == alphabet[k]);
            }
        }

        for p in self.alphabet_size * self.degree as usize..=stop {
            x[p] = x[self.matrices[p].a11 as usize] * x[self.matrices[p].a22 as usize]
                - x[self.matrices[p].a12 as usize] * x[self.matrices[p].a21 as usize];
        }
    }
}

fn tuple_index(letter_counts: &[usize]) -> usize {
    if letter_counts.len() == 2 {
        let s = letter_counts[0] + letter_counts[1];
        return ((s * (s + 1)) >> 1) + letter_counts[1];
    }
    let mut index = 0;
    let mut n = 0;

    for k in 0..letter_counts.len() {
        n += letter_counts[letter_counts.len() - k - 1];
        if n == 0 {
            continue;
        }
        index += usize::try_from(binomial::<i128>(k + n, n - 1))
            .expect("Failed to convert bigint to usize");
    }

    index
}

#[cfg(test)]
mod test {

    use lyndon_rs::lyndon::Sort;
    use num_rational::Ratio;

    use super::*;

    #[test]
    fn test_construction() {
        let basis = LyndonBasis::<char>::new(2, Sort::Lexicographical);
        let bch_series_generator = BchSeriesGenerator::<char>::new(basis, 5);
        assert_eq!(bch_series_generator.max_degree, 5);
        assert_eq!(bch_series_generator.basis.len(), 14);
        let expected_left_factor = [
            0, // A
            1, // B
            0, // AB
            0, // AAB
            2, // ABB
            0, // AAAB
            0, // AABB
            4, // ABBB
            0, // AAAAB
            0, // AAABB
            3, // AABAB
            0, // AABBB
            2, // ABABB
            7, // ABBBB
        ];
        for (term, expected_term) in bch_series_generator
            .left_factor
            .iter()
            .zip(expected_left_factor.iter())
        {
            assert_eq!(term, expected_term);
        }
        let expected_right_factor = [
            0, // A
            0, // B
            1, // AB
            2, // AAB
            1, // ABB
            3, // AAAB
            4, // AABB
            1, // ABBB
            5, // AAAAB
            6, // AAABB
            2, // AABAB
            7, // AABBB
            4, // ABABB
            1, // ABBBB
        ];
        for (term, expected_term) in bch_series_generator
            .right_factor
            .iter()
            .zip(expected_right_factor.iter())
        {
            assert_eq!(term, expected_term);
        }

        let expected_multi_degree_indices = [
            1,  // A
            2,  // B
            4,  // AB
            7,  // AAB
            8,  // ABB
            11, // AAAB
            12, // AABB
            13, // ABBB
            16, // AAAAB
            17, // AAABB
            17, // AABAB
            18, // AABBB
            18, // ABABB
            19, // ABBBB
        ];
        for (term, expected_term) in bch_series_generator
            .multi_degree
            .iter()
            .zip(expected_multi_degree_indices.iter())
        {
            assert_eq!(term, expected_term);
        }
    }

    #[test]
    fn test_goldberg_series() {
        let basis = LyndonBasis::<char>::new(2, Sort::Lexicographical);
        let bch_series_generator = BchSeriesGenerator::<char>::new(basis, 5);
        let basis = basis.generate_basis(5);
        let goldberg_coefficients = bch_series_generator
            .generate_goldberg_coefficient_numerators::<i128>()
            .into_iter()
            .zip(basis.iter())
            .map(|(numerator, word)| {
                let n = word.goldberg().iter().sum();
                let denominator = FACTORIALS[n] * bch_denominator::<i128>(n);
                Ratio::new(numerator, denominator)
            })
            .collect::<Vec<_>>();
        let expected_goldberg_coefficients = vec![
            Ratio::new(1, 1),    // A
            Ratio::new(1, 1),    // B
            Ratio::new(1, 2),    // AB
            Ratio::new(1, 12),   // AAB
            Ratio::new(1, 12),   // ABB
            Ratio::default(),    // AAAB
            Ratio::new(1, 24),   // AABB
            Ratio::default(),    // ABBB
            Ratio::new(-1, 720), // AAAAB
            Ratio::new(1, 180),  // AAABB
            Ratio::new(-1, 120), // AABAB
            Ratio::new(1, 180),  // AABBB
            Ratio::new(-1, 120), // ABABB
            Ratio::new(-1, 720), // ABBBB
        ];
        for (term, expected_term) in goldberg_coefficients
            .iter()
            .zip(expected_goldberg_coefficients.iter())
        {
            assert_eq!(term, expected_term);
        }
    }

    #[test]
    fn test_bch_series() {
        let basis = LyndonBasis::<char>::new(2, Sort::Lexicographical);
        let bch_series_generator = BchSeriesGenerator::<char>::new(basis, 5);
        let bch_coefficients = bch_series_generator.generate_bch_coefficients::<Ratio<i128>>();
        let expected_bch_coefficients = vec![
            Ratio::new(1, 1),
            Ratio::new(1, 1),
            Ratio::new(1, 2),
            Ratio::new(1, 12),
            Ratio::new(1, 12),
            Ratio::default(),
            Ratio::new(1, 24),
            Ratio::default(),
            Ratio::new(-1, 720),
            Ratio::new(1, 180),
            Ratio::new(1, 360),
            Ratio::new(1, 180),
            Ratio::new(1, 120),
            Ratio::new(-1, 720),
        ];

        for (term, expected_term) in bch_coefficients
            .iter()
            .zip(expected_bch_coefficients.iter())
        {
            assert_eq!(term, expected_term);
        }
    }
}
