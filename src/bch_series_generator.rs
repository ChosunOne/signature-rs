use crate::{
    BCHCoefficientGenerator, Int, LieSeriesGenerator,
    bch::{bch_denominator, goldberg_coeff_numerator},
    bernoulli_sequence, binomial,
    constants::FACTORIALS,
    lie_series::LieSeries,
    lyndon::{Generator, LyndonBasis, LyndonWord, Topological},
    rooted_tree::{GraphPartitionTable, RootedTree},
};
use bitvec::prelude::*;
use num_rational::Ratio;
use rayon::prelude::*;
use std::{collections::HashMap, hash::Hash, marker::PhantomData};

#[cfg(feature = "progress")]
mod progress_imports {
    pub use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
}

#[cfg(feature = "progress")]
pub use progress_imports::*;

/// A Lie series of words with alphabet size `N` and generators of type `T`
pub struct BchSeriesGenerator<const N: usize = 2, T: Generator<Letter = T> + Send + Sync = u8> {
    /// Maximum length of the words
    pub max_degree: usize,
    /// The words of the series
    pub basis: Vec<LyndonWord<N, T>>,
    /// The index of the left factor for a given word index `i`
    pub left_factor: Vec<usize>,
    /// The index of the right factor for a given word index `i`
    pub right_factor: Vec<usize>,
    /// The length of a given word `i`
    pub word_lengths: Vec<usize>,
    /// The starting indices of words with degree `n`
    pub index_of_degree: Vec<usize>,
    /// The multi-degree index for a given word index `i`
    pub multi_degree: Vec<usize>,
}

impl<const N: usize, T: Generator<Letter = T> + Send + Sync> BchSeriesGenerator<N, T> {
    pub fn new(basis: Vec<LyndonWord<N, T>>) -> Self {
        #[cfg(feature = "progress")]
        let style = ProgressStyle::with_template(
            "[{eta_precise}] [{bar:35.green/white}] {pos:>2}/{len:2} {msg}",
        )
        .unwrap()
        .progress_chars("=>-");

        let max_degree = basis.last().map(|x| x.len()).unwrap_or(0);
        let number_of_words_per_degree =
            LyndonBasis::<N, T>::number_of_words_per_degree(max_degree);

        let mut left_factor = vec![0; basis.len()];
        for i in 0..N {
            left_factor[i] = i;
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
                    for k in index_of_degree[l.len() - 1]..=index_of_degree[l.len()] - 1 {
                        if basis[k] == l {
                            *lf = k;
                            break;
                        }
                    }
                    for k in index_of_degree[r.len() - 1]..=index_of_degree[r.len()] - 1 {
                        if basis[k] == r {
                            *rf = k;
                            break;
                        }
                    }
                }
                *wl = word.len();
            });

        #[cfg(feature = "progress")]
        pb.finish();

        let alphabet = T::alphabet::<N>();
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
            let mut letter_counts = [0usize; N];
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
            max_degree,
            basis,
            left_factor,
            right_factor,
            word_lengths,
            index_of_degree,
            multi_degree,
        }
    }

    pub fn generate_goldberg_coefficient_numerators<U: Int + Send + Sync>(&self) -> Vec<U> {
        #[cfg(feature = "progress")]
        let style = ProgressStyle::with_template(
            "[{eta_precise}] [{bar:35.green/white}] {pos:>2}/{len:2} {msg}",
        )
        .unwrap()
        .progress_chars("=>-");

        let mut goldberg_coefficient_numerators = vec![U::from(0); self.basis.len()];

        let alphabet = T::alphabet::<N>();
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
                *numer = goldberg_coeff_numerator(goldberg_partition, a_first)
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
}

impl<const N: usize, T: Generator<Letter = T> + Send + Sync> BCHCoefficientGenerator
    for BchSeriesGenerator<N, T>
{
    fn generate_bch_coefficients<U: Int + Send + Sync>(&self) -> Vec<Ratio<U>> {
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
        let alphabet = T::alphabet::<N>();

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

            let mut matrix_tree = MatrixTree::<N, T>::new(self.max_degree as u8);
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

            let mut X = vec![0_i64; matrix_tree.matrices.len()];
            let mut stop = 0;
            for x in 0..words_in_group.len() {
                let i = words_in_group[x];
                stop = if group_matrix_indices[x] > stop {
                    group_matrix_indices[x]
                } else {
                    stop
                };

                let (num_right_factors, right_factors) = self.right_factors(i, self.max_degree);
                let word = &self.basis[i];
                matrix_tree.run(&mut X, word, stop);

                let mut index_of_non_initial_generator = 0;
                while word.letters[index_of_non_initial_generator] == alphabet[0] {
                    index_of_non_initial_generator += 1;
                }

                for y in 0..x {
                    let j = words_in_group[y];
                    let (previous_num_right_factors, previous_right_factors) =
                        self.right_factors(j, self.max_degree);
                    if index_of_non_initial_generator >= previous_num_right_factors {
                        let d = X[group_matrix_indices[y]];
                        if d != 0 {
                            for k in 0..=previous_num_right_factors.min(num_right_factors) {
                                let previous_word_coefficient =
                                    bch_coefficient_numerators[previous_right_factors[k]].clone();
                                // SAFETY: Different mutli-degree groups
                                // are disjoint, so this is safe.
                                unsafe {
                                    let ptr = coeffs_addr as *mut U;
                                    let coeffs = std::slice::from_raw_parts_mut(ptr, len);
                                    coeffs[right_factors[k]] -=
                                        U::from(d) * previous_word_coefficient;
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
                let denominator = U::from(FACTORIALS[n]) * bch_denominator::<U>(n);
                Ratio::new(numerator, denominator)
            })
            .collect()
    }
}

impl<const N: usize, T: Generator<Letter = T> + Send + Sync, U: Hash + Int + Send + Sync>
    LieSeriesGenerator<N, T, U> for BchSeriesGenerator<N, T>
{
    fn generate_lie_series(&self) -> LieSeries<N, T, U> {
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

pub struct MatrixTree<const N: usize, T: Generator<Letter = T>> {
    /// A vector of matrices
    matrices: Vec<Matrix2x2>,
    /// The maximum degree of the words
    degree: u8,
    /// Maps keys to indices in `matrices`.
    memoization_cache: HashMap<usize, usize>,
    _phantom: PhantomData<T>,
}

impl<const N: usize, T: Generator<Letter = T>> MatrixTree<N, T> {
    pub fn new(degree: u8) -> Self {
        let mut matrices = Vec::with_capacity(N * degree as usize);
        let mut memoization_cache = HashMap::new();

        for i in 0..N {
            for j in 0..degree {
                let matrix = Matrix2x2::default();
                let key = ((i as usize) << 8) | j as usize;
                memoization_cache.insert(key, matrices.len());
                matrices.push(matrix);
            }
        }

        Self {
            matrices,
            degree,
            memoization_cache,
            _phantom: PhantomData,
        }
    }

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

    pub fn run(&self, x: &mut [i64], word: &LyndonWord<N, T>, mut stop: usize) {
        let alphabet = T::alphabet::<N>();
        if stop >= self.matrices.len() {
            stop = self.matrices.len() - 1;
        }
        for k in 0..N {
            for i in 0..self.degree {
                x[k * self.degree as usize + i as usize] =
                    if word.letters[i as usize] == alphabet[k] {
                        1
                    } else {
                        0
                    };
            }
        }

        for p in N * self.degree as usize..=stop {
            x[p] = x[self.matrices[p].a11 as usize] * x[self.matrices[p].a22 as usize]
                - x[self.matrices[p].a12 as usize] * x[self.matrices[p].a21 as usize];
        }
    }
}

fn tuple_index<const N: usize>(letter_counts: &[usize; N]) -> usize {
    if N == 2 {
        let s = letter_counts[0] as usize + letter_counts[1] as usize;
        return (((s * (s + 1)) >> 1) + letter_counts[1] as usize) as usize;
    }
    let mut index = 0;
    let mut n = 0;

    for k in 0..N {
        n += letter_counts[N - k - 1];
        if n == 0 {
            continue;
        }
        index += usize::try_from(binomial::<i128>(k + n as usize, (n - 1) as usize).clone())
            .expect("Failed to convert bigint to usize");
    }

    index
}

pub struct RootedTreeLieSeries<T: Generator<Letter = T> = u8, U: Int = i128> {
    graph_partition_table: GraphPartitionTable<T>,
    is_computed_z: BitVec,
    z: Vec<Ratio<U>>,
    bernoulli: Vec<Ratio<U>>,
    x_minus_y: Vec<Ratio<U>>,
    x_plus_y: Vec<Ratio<U>>,
    prime: Vec<usize>,
    d_prime: Vec<usize>,
    sigma: Vec<usize>,
    kappa: Vec<usize>,
    m_n: usize,
    adjoint_cache: HashMap<(usize, usize), Ratio<U>>,
    progress: u64,
}

impl<T: Generator<Letter = T>, U: Int> RootedTreeLieSeries<T, U> {
    pub fn new(n: usize) -> Self {
        let t_n = LyndonBasis::<2, T, Topological>::generate_basis(n)
            .into_iter()
            .map(|w| RootedTree::from(w))
            .collect::<Vec<_>>();
        let m_n = t_n.len();
        let graph_partition_table = GraphPartitionTable::new(t_n);
        let is_computed_z = bitvec![0; graph_partition_table.tm_n()];
        let bernoulli = bernoulli_sequence(n);
        let z = vec![Ratio::default(); graph_partition_table.tm_n()];
        let mut x_minus_y = vec![Ratio::default(); graph_partition_table.tm_n()];
        x_minus_y[0] = Ratio::new(1.into(), 1.into());
        x_minus_y[1] = Ratio::new((-1).into(), 1.into());
        let mut x_plus_y = vec![Ratio::default(); graph_partition_table.tm_n()];
        x_plus_y[0] = Ratio::new(1.into(), 1.into());
        x_plus_y[1] = Ratio::new(1.into(), 1.into());
        let prime = vec![0; m_n];
        let d_prime = vec![0; m_n];
        let mut sigma = vec![0_usize; m_n];
        sigma[0] = 1;
        sigma[1] = 1;
        let mut kappa = vec![0; m_n];
        kappa[0] = 1;
        kappa[1] = 1;
        let adjoint_cache = HashMap::new();

        Self {
            graph_partition_table,
            is_computed_z,
            bernoulli,
            z,
            x_minus_y,
            x_plus_y,
            prime,
            d_prime,
            sigma,
            kappa,
            m_n,
            adjoint_cache,
            progress: 0,
        }
    }

    fn lie_bracket(&self, alpha: &[Ratio<U>], beta: &[Ratio<U>], i: usize) -> Ratio<U> {
        let mut sum = Ratio::<U>::default();
        for &(j, k) in self.graph_partition_table.partitions(i).partitions.iter() {
            sum += &alpha[j] * &beta[k] - &alpha[k] * &beta[j];
        }

        sum
    }

    fn adjoint_operator(&mut self, i: usize, power: usize) -> Ratio<U> {
        if power == 0 {
            return self.x_plus_y[i].clone();
        }
        if let Some(v) = self.adjoint_cache.get(&(i, power)) {
            return v.clone();
        }
        if power == 1 {
            let result = self.lie_bracket(&self.z, &self.x_plus_y, i);
            self.adjoint_cache.insert((i, power), result.clone());
            return result;
        }
        let mut sum = Ratio::default();
        let partitions = self
            .graph_partition_table
            .partitions(i)
            .partitions
            .iter()
            .copied()
            .collect::<Vec<_>>();
        for (j, k) in partitions {
            sum += self.z[j].clone() * self.adjoint_operator(k, power - 1)
                - self.z[k].clone() * self.adjoint_operator(j, power - 1);
        }
        self.adjoint_cache.insert((i, power), sum.clone());

        sum
    }

    fn compute_z(&mut self, i: usize) -> Ratio<U> {
        if i == 0 || i == 1 {
            self.is_computed_z.set(i, true);
            return Ratio::new(1.into(), 1.into());
        }
        // 1/2 [X-Y, Z] (u_i)
        let s_ui = self.graph_partition_table.partitions(i).clone();
        for &(j, k) in &s_ui.partitions {
            // Check if there are any Z_u terms that we need to complete first
            if !self.is_computed_z[j] {
                self.z[j] = self.compute_z(j);
                self.progress += 1;
            }
            if !self.is_computed_z[k] {
                self.z[k] = self.compute_z(k);
                self.progress += 1;
            }
        }
        let left_term =
            Ratio::new(1.into(), 2.into()) * self.lie_bracket(&self.x_minus_y, &self.z, i);

        // Bernoulli term
        let mut bernoulli_term = Ratio::default();
        let degree = self.graph_partition_table.degree(i);
        for p in 1..=((degree - 1) / 2) {
            let bernoulli_coef =
                &self.bernoulli[2 * p] * Ratio::new(1.into(), (FACTORIALS[2 * p]).into());
            let adjoint_term = self.adjoint_operator(i, 2 * p);
            bernoulli_term += bernoulli_coef * adjoint_term;
        }
        self.is_computed_z.set(i, true);
        let result = (left_term + bernoulli_term)
            * Ratio::new(
                1.into(),
                (self.graph_partition_table.degree(i) as i8).into(),
            );
        result
    }

    /// Generates the truncated Baker-Campbell-Hausdorff series of a Lyndon basis
    pub fn generate_coefficients(&mut self) -> Vec<Ratio<U>> {
        #[cfg(feature = "progress")]
        let pb = ProgressBar::new(self.graph_partition_table.tm_n() as u64);
        #[cfg(feature = "progress")]
        {
            pb.set_message("Generating BCH Coefficients ");
            pb.tick();
        }
        for i in 0..self.m_n {
            if !self.is_computed_z[i] {
                self.z[i] = self.compute_z(i);
                self.progress += 1;
                if i > 1 {
                    self.prime[i] = self.graph_partition_table.partitions(i)[0].0;
                    self.d_prime[i] = self.graph_partition_table.partitions(i)[0].1;
                    self.kappa[i] = {
                        if self.d_prime[self.prime[i]] != self.d_prime[i] {
                            1
                        } else {
                            self.kappa[self.prime[i]] + 1
                        }
                    };
                    self.sigma[i] =
                        self.kappa[i] * self.sigma[self.prime[i]] * self.sigma[self.d_prime[i]];
                }
            }

            #[cfg(feature = "progress")]
            pb.set_position(self.progress);
        }
        #[cfg(feature = "progress")]
        pb.finish();

        self.z[..self.m_n]
            .iter()
            .zip(self.sigma.iter())
            .map(|(z, &sig)| z * Ratio::new(1.into(), (sig as u64).into()))
            .collect()
    }
}

#[cfg(test)]
mod test {
    use crate::lyndon::Lexicographical;

    use super::*;

    #[test]
    fn test_lie_series_construction() {
        let basis = LyndonBasis::<2, char, Lexicographical>::generate_basis(5);
        let lie_series = BchSeriesGenerator::<2, char>::new(basis);
        assert_eq!(lie_series.max_degree, 5);
        assert_eq!(lie_series.basis.len(), 14);
        let expected_left_factor = vec![
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
        for (i, (term, expected_term)) in lie_series
            .left_factor
            .iter()
            .zip(expected_left_factor.iter())
            .enumerate()
        {
            dbg!(i);
            assert_eq!(term, expected_term);
        }
        let expected_right_factor = vec![
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
        dbg!(&lie_series.right_factor);
        for (i, (term, expected_term)) in lie_series
            .right_factor
            .iter()
            .zip(expected_right_factor.iter())
            .enumerate()
        {
            dbg!(i);
            assert_eq!(term, expected_term);
        }

        let expected_multi_degree_indices = vec![
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
        for (i, (term, expected_term)) in lie_series
            .multi_degree
            .iter()
            .zip(expected_multi_degree_indices.iter())
            .enumerate()
        {
            dbg!(i);
            assert_eq!(term, expected_term);
        }
    }

    #[test]
    fn test_lie_series_goldberg_series() {
        let basis = LyndonBasis::<2, char, Lexicographical>::generate_basis(5);
        let lie_series = BchSeriesGenerator::<2, char>::new(basis.clone());
        let goldberg_coefficients = lie_series
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
        dbg!(&goldberg_coefficients);
        for (i, (term, expected_term)) in goldberg_coefficients
            .iter()
            .zip(expected_goldberg_coefficients.iter())
            .enumerate()
        {
            dbg!(i);
            assert_eq!(term, expected_term);
        }
    }

    #[test]
    fn test_lie_series_bch_series() {
        let basis = LyndonBasis::<2, char, Lexicographical>::generate_basis(5);
        let lie_series = BchSeriesGenerator::<2, char>::new(basis);
        let bch_coefficients = lie_series.generate_bch_coefficients();
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

        dbg!(&bch_coefficients);
        for (i, (term, expected_term)) in bch_coefficients
            .iter()
            .zip(expected_bch_coefficients.iter())
            .enumerate()
        {
            dbg!(i);
            assert_eq!(term, expected_term)
        }
    }

    #[test]
    fn test_rooted_tree_bch_series() {
        let mut generator = RootedTreeLieSeries::<char>::new(5);
        let series = generator.generate_coefficients();
        let expected_z_ui_series = vec![
            Ratio::new(1.into(), 1.into()),
            Ratio::new(1.into(), 1.into()),
            Ratio::new(1.into(), 2.into()),
            Ratio::new(1.into(), 12.into()),
            Ratio::new(1.into(), 12.into()),
            Ratio::default(),
            Ratio::new(1.into(), 24.into()),
            Ratio::default(),
            Ratio::new((-1).into(), 720.into()),
            Ratio::new(1.into(), 360.into()),
            Ratio::new(1.into(), 120.into()),
            Ratio::new(1.into(), 180.into()),
            Ratio::new(1.into(), 180.into()),
            Ratio::new((-1).into(), 720.into()),
        ];
        for (term, expected_term) in series.iter().zip(&expected_z_ui_series) {
            assert_eq!(term, expected_term);
        }
    }
}
