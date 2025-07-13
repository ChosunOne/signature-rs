use std::{collections::HashMap, marker::PhantomData};

use crate::{
    bch::goldberg_coeff,
    binomial,
    lyndon::{Generator, LyndonBasis, LyndonWord},
};

#[cfg(feature = "progress")]
use indicatif::ProgressBar;
#[cfg(feature = "progress")]
use indicatif::{MultiProgress, ProgressStyle};
use num_integer::Integer;
use num_rational::Ratio;

/// A Lie series of words with alphabet size `N` and generators of type `T`
pub struct LieSeries<const N: usize, T: Generator<Letter = T>, U: Clone + Integer = i128> {
    /// Maximum length of the words
    pub max_degree: usize,
    /// The words of the series
    pub words: Vec<LyndonWord<N, T>>,
    /// The index of the left factor for a given word index `i`
    pub left_factor: Vec<usize>,
    /// The index of the right factor for a given word index `i`
    pub right_factor: Vec<usize>,
    /// The length of a given word `i`
    pub word_lengths: Vec<usize>,
    /// The starting indices of words with degree `n`
    pub index_of_degree: Vec<usize>,
    /// The Goldberg coefficients for a given word index `i`
    pub goldberg_coefficients: Vec<Ratio<U>>,
    /// The Baker-Campbell-Hausdorff coefficients for a given word index `i`
    pub bch_coefficients: Vec<Ratio<U>>,
    /// The multi-degree index for a given word index `i`
    pub multi_degree: Vec<usize>,
}

impl<const N: usize, T: Generator<Letter = T>, U: Clone + Integer> LieSeries<N, T, U> {
    pub fn new(max_degree: usize) -> Self {
        #[cfg(feature = "progress")]
        let style = ProgressStyle::with_template(
            "[{eta_precise}] {bar:35.green/white} {pos:>2}/{len:2} {msg}",
        )
        .unwrap()
        .progress_chars("=>-");
        let number_of_words_per_degree =
            LyndonBasis::<N, T>::number_of_words_per_degree(max_degree);
        let words = LyndonBasis::<N, T>::generate_basis(max_degree, false);
        let mut word_index_map = HashMap::new();
        for (i, word) in words.iter().enumerate() {
            word_index_map.insert(word, i);
        }

        let mut left_factor = vec![0; words.len()];
        for i in 0..N {
            left_factor[i] = i;
        }
        let mut right_factor = vec![0; words.len()];
        let mut word_lengths = vec![0; words.len()];
        let mut index_of_degree = vec![0; max_degree + 1];

        for n in 1..=max_degree {
            index_of_degree[n] = index_of_degree[n - 1] + number_of_words_per_degree[n - 1];
        }

        #[cfg(feature = "progress")]
        let pb = ProgressBar::new(words.len() as u64).with_style(style.clone());
        #[cfg(feature = "progress")]
        {
            pb.set_message("Factorizing Lyndon words ");
            pb.tick();
        }
        for (i, word) in words.iter().enumerate() {
            if word.len() > 1 {
                let (l, r) = word.factorize();
                for k in index_of_degree[l.len() - 1]..=index_of_degree[l.len()] - 1 {
                    if words[k] == l {
                        left_factor[i] = k;
                        break;
                    }
                }
                for k in index_of_degree[r.len() - 1]..=index_of_degree[r.len()] - 1 {
                    if words[k] == r {
                        right_factor[i] = k;
                        break;
                    }
                }
            }
            word_lengths[i] = word.len();
            #[cfg(feature = "progress")]
            pb.inc(1);
        }

        #[cfg(feature = "progress")]
        pb.finish();

        let mut goldberg_coefficients = vec![Ratio::default(); words.len()];

        let alphabet = T::alphabet::<N>();
        #[cfg(feature = "progress")]
        {
            pb.set_length(0);
            pb.set_message("Calculating Goldberg Coefficients ");
            pb.tick();
        }
        for (i, word) in words.iter().enumerate() {
            let goldberg_partition = word.goldberg();
            let a_first = word.letters[0] == alphabet[0];
            goldberg_coefficients[i] = goldberg_coeff(goldberg_partition, a_first);
            #[cfg(feature = "progress")]
            pb.inc(1);
        }

        #[cfg(feature = "progress")]
        pb.finish();

        let mut bch_coefficients = goldberg_coefficients.clone();

        let mut multi_degree = vec![0; words.len()];
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
        for i in 0..words.len() {
            let mut letter_counts = [0usize; N];
            let word = &words[i];
            for j in 0..word_lengths[i] {
                letter_counts[alphabet_index[&word.letters[j]]] += 1;
            }
            multi_degree[i] = tuple_index(&letter_counts);
            #[cfg(feature = "progress")]
            pb.inc(1);
        }
        #[cfg(feature = "progress")]
        {
            pb.finish();
            drop(pb);
        }

        // Loop over the words of max degree
        let start = index_of_degree[max_degree - 1];
        let end = index_of_degree[max_degree] - 1;
        let h1 = multi_degree[start];
        let h2 = multi_degree[end];

        #[cfg(feature = "progress")]
        let mb = MultiProgress::new();
        #[cfg(feature = "progress")]
        mb.println("Calculating BCH Coefficients").unwrap();
        #[cfg(feature = "progress")]
        let p1 = mb.add(ProgressBar::new((h2 - h1 + 1) as u64).with_style(style.clone()));
        #[cfg(feature = "progress")]
        let p2 = mb.add(ProgressBar::new(0).with_style(style.clone()));
        #[cfg(feature = "progress")]
        let p3 = mb.add(ProgressBar::new(0).with_style(style.clone()));
        #[cfg(feature = "progress")]
        let p4 = mb.add(ProgressBar::new(0).with_style(style.clone()));

        #[cfg(feature = "progress")]
        p1.set_message("Calculating by multi-degree ");

        for h in h1..=h2 {
            let words_in_group = words[start..=end]
                .iter()
                .enumerate()
                .filter(|&(i, _)| multi_degree[i + start] == h)
                .map(|(i, _)| i + start)
                .collect::<Vec<_>>();

            let mut matrix_tree = MatrixTree::<N, T>::new(max_degree as u8);
            let mut group_matrix_indices = Vec::with_capacity(words_in_group.len());
            #[cfg(feature = "progress")]
            {
                p2.set_position(0);
                p2.set_length(words_in_group.len() as u64);
                p2.set_message("Finding words in group ");
                p2.tick();
            }

            for &j in &words_in_group {
                let (num_factors, right_factors) =
                    get_right_factors(j, max_degree, &left_factor, &right_factor);

                let rightmost_factor = right_factors[num_factors];

                group_matrix_indices.push(matrix_tree.append(
                    rightmost_factor,
                    num_factors as u8,
                    &left_factor,
                    &right_factor,
                    &word_lengths,
                ));
                #[cfg(feature = "progress")]
                p2.inc(1);
            }

            #[cfg(feature = "progress")]
            {
                p2.set_position(0);
                p2.set_length(words_in_group.len() as u64);
                p2.set_message("Matrix Tree Runs ");
            }
            let mut X = vec![0_isize; matrix_tree.matrices.len()];
            let mut stop = 0;
            for x in 0..words_in_group.len() {
                let i = words_in_group[x];
                stop = if group_matrix_indices[x] > stop {
                    group_matrix_indices[x]
                } else {
                    stop
                };

                let (num_right_factors, right_factors) =
                    get_right_factors(i, max_degree, &left_factor, &right_factor);
                let word = &words[i];
                matrix_tree.run(&mut X, word, stop);

                let mut index_of_non_initial_generator = 0;
                while word.letters[index_of_non_initial_generator] == alphabet[0] {
                    index_of_non_initial_generator += 1;
                }

                #[cfg(feature = "progress")]
                {
                    p3.set_position(0);
                    p3.set_length(x as u64);
                    p3.set_message("Iterating words in group ");
                    p3.tick();
                }
                for y in 0..x {
                    let j = words_in_group[y];
                    let (previous_num_right_factors, previous_right_factors) =
                        get_right_factors(j, max_degree, &left_factor, &right_factor);
                    if index_of_non_initial_generator >= previous_num_right_factors {
                        let d = X[group_matrix_indices[y]];
                        if d != 0 {
                            #[cfg(feature = "progress")]
                            {
                                p4.set_position(0);
                                p4.set_length(
                                    previous_num_right_factors.min(num_right_factors) as u64
                                );
                                p4.set_message("Generating BCH terms");
                                p4.tick();
                            }

                            for k in 0..=previous_num_right_factors.min(num_right_factors) {
                                let previous_word_coefficient =
                                    bch_coefficients[previous_right_factors[k]].clone();
                                bch_coefficients[right_factors[k]] -=
                                    Ratio::new(d.try_into().unwrap(), 1.into())
                                        * previous_word_coefficient;
                                #[cfg(feature = "progress")]
                                {
                                    p4.inc(1);
                                    p3.tick();
                                    p2.tick();
                                    p1.tick();
                                }
                            }
                        }
                    }
                    #[cfg(feature = "progress")]
                    p3.inc(1);
                }
                #[cfg(feature = "progress")]
                p2.inc(1);
            }
            #[cfg(feature = "progress")]
            p1.inc(1);
        }

        #[cfg(feature = "progress")]
        {
            p1.finish();
            p2.finish();
            p3.finish();
            p4.finish();
        }

        Self {
            max_degree,
            words,
            left_factor,
            right_factor,
            word_lengths,
            index_of_degree,
            goldberg_coefficients,
            bch_coefficients,
            multi_degree,
        }
    }
}

// TODO: Make this a method of the lie series
fn get_right_factors(
    i: usize,
    kmax: usize,
    left_factors: &[usize],
    right_factors: &[usize],
) -> (usize, Vec<usize>) {
    let mut j = vec![0; kmax];
    let mut k = 0;
    j[0] = i;
    let mut l = i;
    while k < kmax && left_factors[l] == 0 {
        k += 1;
        l = right_factors[l];
        j[k] = l;
    }
    (k, j)
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

    pub fn run(&self, x: &mut [isize], word: &LyndonWord<N, T>, mut stop: usize) {
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
        #[cfg(not(feature = "bigint"))]
        {
            index += *binomial(k + n as usize, (n - 1) as usize).numer() as usize;
        }
        #[cfg(feature = "bigint")]
        {
            index += usize::try_from(binomial(k + n as usize, (n - 1) as usize).numer().clone())
                .expect("Failed to convert bigint to usize");
        }
    }

    index
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_lie_series_construction() {
        let lie_series = LieSeries::<2, char>::new(5);
        assert_eq!(lie_series.max_degree, 5);
        assert_eq!(lie_series.words.len(), 14);
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

        let expected_goldberg_coefficients = vec![
            Rational::new(1, 1),    // A
            Rational::new(1, 1),    // B
            Rational::new(1, 2),    // AB
            Rational::new(1, 12),   // AAB
            Rational::new(1, 12),   // ABB
            Rational::default(),    // AAAB
            Rational::new(1, 24),   // AABB
            Rational::default(),    // ABBB
            Rational::new(-1, 720), // AAAAB
            Rational::new(1, 180),  // AAABB
            Rational::new(-1, 120), // AABAB
            Rational::new(1, 180),  // AABBB
            Rational::new(-1, 120), // ABABB
            Rational::new(-1, 720), // ABBBB
        ];
        dbg!(&lie_series.goldberg_coefficients);
        for (i, (term, expected_term)) in lie_series
            .goldberg_coefficients
            .iter()
            .zip(expected_goldberg_coefficients.iter())
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

        let expected_bch_coefficients = vec![
            Rational::new(1, 1),
            Rational::new(1, 1),
            Rational::new(1, 2),
            Rational::new(1, 12),
            Rational::new(1, 12),
            Rational::default(),
            Rational::new(1, 24),
            Rational::default(),
            Rational::new(-1, 720),
            Rational::new(1, 180),
            Rational::new(1, 360),
            Rational::new(1, 180),
            Rational::new(1, 120),
            Rational::new(-1, 720),
        ];

        for (i, (term, expected_term)) in lie_series
            .bch_coefficients
            .iter()
            .zip(expected_bch_coefficients.iter())
            .enumerate()
        {
            dbg!(i);
            assert_eq!(term, expected_term)
        }
    }
}
