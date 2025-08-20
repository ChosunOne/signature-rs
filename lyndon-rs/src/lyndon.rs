use std::{
    collections::HashMap,
    fmt::{Debug, Display},
    hash::Hash,
    marker::PhantomData,
    ops::Mul,
    str::FromStr,
};

use thiserror::Error;

use crate::generators::Generator;

/// Sorting method for ordering Lyndon words in a basis.
#[derive(Copy, Clone, Default, Debug, PartialEq, Eq)]
pub enum Sort {
    /// Standard dictionary ordering based on lexicographic comparison.
    #[default]
    Lexicographical,
    /// Ordering based on the structure of Lyndon word factorizations.
    Topological,
}

/// A basis of Lyndon words over an alphabet of type `T`.
///
/// A Lyndon basis provides an ordered set of Lyndon words that can be used
/// as a basis for the free Lie algebra over a given alphabet.
pub struct LyndonBasis<T> {
    /// The size of the alphabet used to generate Lyndon words.
    pub alphabet_size: usize,
    /// The sorting method used to order words in the basis.
    sort: Sort,
    /// Phantom data to associate the basis with generator type `T`.
    _generator: PhantomData<T>,
}

impl<T> Copy for LyndonBasis<T> {}

impl<T> Clone for LyndonBasis<T> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<T> Default for LyndonBasis<T> {
    fn default() -> Self {
        Self {
            alphabet_size: usize::default(),
            sort: Sort::default(),
            _generator: PhantomData,
        }
    }
}

impl<T> Debug for LyndonBasis<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("LyndonBasis")
            .field("alphabet_size", &self.alphabet_size)
            .field("sort", &self.sort)
            .finish()
    }
}

impl<T> PartialEq for LyndonBasis<T> {
    fn eq(&self, other: &Self) -> bool {
        self.alphabet_size == other.alphabet_size
            && self.sort == other.sort
            && self._generator == other._generator
    }
}

impl<T> Eq for LyndonBasis<T> {}

impl<T> LyndonBasis<T> {
    /// Creates a new Lyndon basis with the specified alphabet size and sorting method.
    #[must_use]
    pub fn new(alphabet_size: usize, sort: Sort) -> Self {
        Self {
            alphabet_size,
            sort,
            _generator: PhantomData,
        }
    }

    /// Calculates the number of Lyndon words of each degree up to `max_degree`.
    ///
    /// Uses the Möbius function to compute the exact count of Lyndon words
    /// for each length from 1 to `max_degree`.
    #[must_use]
    pub fn number_of_words_per_degree(&self, max_degree: usize) -> Vec<usize> {
        let mu = moebius_mu(max_degree);
        let mut words_per_degree = vec![0; max_degree];
        for n in 1..=max_degree {
            let mut d = 1_i64;
            let mut h = 0_i64;
            while d * d < n as i64 {
                let quot = n as i64 / d;
                let rem = n as i64 % d;
                if rem == 0 {
                    h += mu[(d - 1) as usize] * (self.alphabet_size as i64).pow(quot as u32)
                        + mu[(quot - 1) as usize] * (self.alphabet_size as i64).pow(d as u32);
                }
                d += 1;
            }
            if d * d == n as i64 {
                h += mu[(d - 1) as usize] * (self.alphabet_size as i64).pow(d as u32);
            }
            words_per_degree[n - 1] = h as usize / n;
        }

        words_per_degree
    }
}

impl<T: Generator + Clone + Eq + Hash + Ord> LyndonBasis<T> {
    /// Generates all Lyndon words up to the specified maximum length.
    ///
    /// Returns a vector of Lyndon words ordered according to the basis sorting method.
    #[must_use]
    pub fn generate_basis(&self, max_length: usize) -> Vec<LyndonWord<T>> {
        let total_words: usize = self.number_of_words_per_degree(max_length).iter().sum();

        let alphabet = T::alphabet(self.alphabet_size);
        let letter_index: HashMap<_, _> = alphabet
            .iter()
            .cloned()
            .enumerate()
            .map(|(i, v)| (v, i))
            .collect();
        let mut basis = Vec::with_capacity(total_words);
        if max_length == 0 {
            return basis;
        }
        let mut w: Vec<T> = vec![];

        loop {
            if let Some(l) = w.last_mut() {
                let i = letter_index[&*l];
                *l = alphabet[i + 1].clone();
            } else {
                w.push(alphabet[0].clone());
            }

            if let Some(l) = w.last()
                && w.len() <= max_length
                && l.clone() <= alphabet[self.alphabet_size - 1]
            {
                basis.push(LyndonWord::try_from(w.clone()).expect("To make a lyndon word"));

                let m = w.len();
                while w.len() < max_length {
                    w.push(w[w.len() % m].clone());
                }
            }

            while !w.is_empty() && *w.last().unwrap() >= alphabet[self.alphabet_size - 1] {
                w.pop();
            }

            if w.is_empty() {
                break;
            }
        }
        basis.sort_by_key(|word| (word.len(), word.letters.clone()));

        match self.sort {
            Sort::Lexicographical => basis,
            Sort::Topological => Self::topological_sort(&basis),
        }
    }

    fn topological_sort(basis: &[LyndonWord<T>]) -> Vec<LyndonWord<T>> {
        let max_length = basis.last().map_or(0, LyndonWord::len);
        let mut sorted_basis = Vec::with_capacity(basis.len());
        let mut sorted_basis_index = HashMap::new();

        for level in 1..=max_length {
            let mut unsorted_words_by_level = basis
                .iter()
                .filter(|word| word.len() == level)
                .collect::<Vec<_>>();

            if level == 1 {
                // Sort lexicographically
                unsorted_words_by_level.sort_by_key(|&w| w);
                for word in unsorted_words_by_level {
                    sorted_basis_index.insert(word.clone(), sorted_basis.len());
                    sorted_basis.push(word.clone());
                }

                continue;
            }

            // Topological Sort
            unsorted_words_by_level.sort_by_key(|w| {
                let (_, i_dp) = w.factorize();
                sorted_basis_index[&i_dp]
            });
            for word in unsorted_words_by_level {
                sorted_basis_index.insert(word.clone(), sorted_basis.len());
                sorted_basis.push(word.clone());
            }
        }

        sorted_basis
    }
}

/// Computes the Möbius function values up to `max_degree`.
///
/// The Möbius function μ(n) is used in the computation of Lyndon word counts
/// and has the following properties:
/// - μ(n) = 1 if n is a square-free positive integer with an even number of prime factors
/// - μ(n) = -1 if n is a square-free positive integer with an odd number of prime factors  
/// - μ(n) = 0 if n has a squared prime factor
#[must_use]
pub fn moebius_mu(max_degree: usize) -> Vec<i64> {
    let mut mu = vec![0; max_degree];
    mu[0] = 1;
    for n in 1..=max_degree / 2 {
        let mu_n = mu[n - 1];
        for i in (2 * n - 1..max_degree).step_by(n) {
            mu[i] -= mu_n;
        }
    }

    mu
}

/// A Lyndon word over an alphabet of type `T`.
///
/// A Lyndon word is a string that is strictly lexicographically smaller than
/// all of its nontrivial rotations. Lyndon words form a basis for the free
/// Lie algebra.
pub struct LyndonWord<T> {
    /// The sequence of letters that make up the Lyndon word.
    pub letters: Vec<T>,
}

impl<T: Display> Display for LyndonWord<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for letter in &self.letters {
            write!(f, "{letter}")?;
        }
        Ok(())
    }
}

impl<T: Debug> Debug for LyndonWord<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("LyndonWord")
            .field("letters", &self.letters)
            .finish()
    }
}

impl<T: Clone> Clone for LyndonWord<T> {
    fn clone(&self) -> Self {
        Self {
            letters: self.letters.clone(),
        }
    }
}

impl<T: PartialEq> PartialEq for LyndonWord<T> {
    fn eq(&self, other: &Self) -> bool {
        self.letters == other.letters
    }
}

impl<T: Eq> Eq for LyndonWord<T> {}

impl<T: PartialEq + PartialOrd> PartialOrd for LyndonWord<T> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.letters.partial_cmp(&other.letters)
    }
}

impl<T: Ord> Ord for LyndonWord<T> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.letters.cmp(&other.letters)
    }
}

impl<T: Hash> Hash for LyndonWord<T> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.letters.hash(state);
    }
}

impl<T> LyndonWord<T> {
    /// Returns the length of the Lyndon word (number of letters).
    #[must_use]
    pub fn len(&self) -> usize {
        self.letters.len()
    }

    /// Returns `true` if the Lyndon word contains no letters.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.letters.is_empty()
    }
}
impl<T: PartialEq> LyndonWord<T> {
    /// Computes the canonical Goldberg representation of the Lyndon word.
    ///
    /// The Goldberg representation expresses a Lyndon word as a partition
    /// based on runs of identical consecutive letters.
    #[must_use]
    pub fn goldberg(&self) -> Vec<usize> {
        if self.letters.is_empty() {
            return vec![];
        }
        let n = self.letters.len();
        let mut run_counts = vec![0; n];
        let mut current_run = 1;

        for i in 1..n {
            if self.letters[i] == self.letters[i - 1] {
                current_run += 1;
            } else {
                run_counts[current_run - 1] += 1;
                current_run = 1;
            }
        }
        run_counts[current_run - 1] += 1;

        let mut partition = Vec::new();

        for i in (0..n).rev() {
            for _ in 0..run_counts[i] {
                partition.push(i + 1);
            }
        }

        partition
    }
}

impl<T: Ord> LyndonWord<T> {
    fn is_lyndon(word: &[T]) -> bool {
        let n = word.len();
        if n == 0 {
            return false;
        }

        let mut i = 0;
        let mut j = 1;
        while j < n {
            match word[i].cmp(&word[j]) {
                std::cmp::Ordering::Equal => {
                    i += 1;
                    j += 1;
                }
                std::cmp::Ordering::Less => {
                    i = 0;
                    j += 1;
                }
                std::cmp::Ordering::Greater => {
                    return false;
                }
            }
        }

        let period = j - i;
        period == n
    }
}

impl<T: Clone + Generator + Ord> LyndonWord<T> {
    /// Computes all right factors of this Lyndon word.
    ///
    /// Right factors are suffixes of the Lyndon word that are themselves Lyndon words.
    #[must_use]
    pub fn right_factors(&self) -> Vec<Self> {
        let alphabet = T::alphabet(1);
        let mut factors = vec![];
        factors.push(self.clone());
        if self.len() > 1 {
            let (v, w) = self.factorize();
            if v.letters[0] == alphabet[0] && v.len() == 1 {
                let w_right_factors = w.right_factors();
                for factor in w_right_factors {
                    factors.push(factor);
                }
            }
        }
        factors
    }
}

impl<T: Clone + Ord> LyndonWord<T> {
    /// Factorizes the Lyndon word into its canonical factorization (v, w).
    ///
    /// Every Lyndon word of length > 1 can be uniquely written as a concatenation
    /// vw where v is a Lyndon word and w is the lexicographically largest Lyndon
    /// suffix of the original word.
    #[must_use]
    pub fn factorize(&self) -> (LyndonWord<T>, LyndonWord<T>) {
        let n = self.letters.len();
        assert!(n > 1, "Word length must be greater than 1.");

        let mut v = vec![];
        let mut w = vec![];
        for split in 1..n {
            if Self::is_lyndon(&self.letters[split..]) {
                v = self.letters[..split].to_vec();
                w = self.letters[split..].to_vec();
                break;
            }
        }
        (
            Self::try_from(v).expect("A factorized lyndon word should produce a lyndon word"),
            Self::try_from(w).expect("A factorized lyndon word should produce a lyndon word"),
        )
    }
}

impl<T: Clone + Ord> Mul for LyndonWord<T> {
    type Output = Result<Self, LyndonWordError>;

    fn mul(self, rhs: Self) -> Self::Output {
        Self::try_from([self.letters, rhs.letters].concat())
    }
}

/// Errors that can occur when working with Lyndon words.
#[derive(Error, Debug)]
pub enum LyndonWordError {
    /// The provided word does not satisfy the Lyndon word property.
    #[error("Word is not a valid Lyndon word")]
    InvalidWord,
    /// The word contains letters that are not part of the expected alphabet.
    #[error("Word uses letters not in the alphabet")]
    InvalidLetter,
}

impl<T: Ord> TryFrom<Vec<T>> for LyndonWord<T> {
    type Error = LyndonWordError;

    fn try_from(value: Vec<T>) -> Result<Self, Self::Error> {
        if !Self::is_lyndon(&value) {
            return Err(LyndonWordError::InvalidWord);
        }
        Ok(Self { letters: value })
    }
}

impl FromStr for LyndonWord<char> {
    type Err = LyndonWordError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::try_from(s.chars().collect::<Vec<_>>())
    }
}

#[cfg(test)]
mod test {

    use rstest::rstest;

    use super::*;

    #[test]
    fn test_make_a_lyndon_word() -> Result<(), LyndonWordError> {
        let letters = "AAB";
        let word = letters.parse::<LyndonWord<char>>()?;
        assert_eq!(&format!("{word}"), "AAB");
        Ok(())
    }

    #[test]
    fn test_factorize_a_lyndon_word() -> Result<(), LyndonWordError> {
        let letters = "AAAAAB";
        let word = letters.parse::<LyndonWord<char>>()?;
        let (v, w) = word.factorize();
        assert_eq!(&format!("{v}"), "A");
        assert_eq!(&format!("{w}"), "AAAAB");
        Ok(())
    }

    #[test]
    fn test_graft_two_lyndon_words() -> Result<(), LyndonWordError> {
        let letters = "AAB";
        let a = letters.parse::<LyndonWord<char>>()?;
        let letters = "AB";
        let b = letters.parse::<LyndonWord<char>>()?;
        let c = (a * b)?;
        assert_eq!(&format!("{c}"), "AABAB");
        Ok(())
    }

    #[test]
    fn test_generate_a_lyndon_basis_1() {
        let basis = LyndonBasis::new(5, Sort::Lexicographical).generate_basis(5);
        assert_eq!(basis.len(), 829);
        for word in &basis {
            assert!(LyndonWord::<u8>::is_lyndon(&word.letters));
        }
    }

    #[test]
    fn test_generate_a_lyndon_basis_2() -> Result<(), LyndonWordError> {
        let basis = LyndonBasis::new(2, Sort::Topological).generate_basis(5);
        let expected_basis = vec![
            LyndonWord::try_from(vec![0])?,
            LyndonWord::try_from(vec![1])?,
            LyndonWord::try_from(vec![0, 1])?,
            LyndonWord::try_from(vec![0, 1, 1])?,
            LyndonWord::try_from(vec![0, 0, 1])?,
            LyndonWord::try_from(vec![0, 1, 1, 1])?,
            LyndonWord::try_from(vec![0, 0, 1, 1])?,
            LyndonWord::try_from(vec![0, 0, 0, 1])?,
            LyndonWord::try_from(vec![0, 1, 1, 1, 1])?,
            LyndonWord::try_from(vec![0, 0, 1, 0, 1])?,
            LyndonWord::try_from(vec![0, 1, 0, 1, 1])?,
            LyndonWord::try_from(vec![0, 0, 1, 1, 1])?,
            LyndonWord::try_from(vec![0, 0, 0, 1, 1])?,
            LyndonWord::try_from(vec![0, 0, 0, 0, 1])?,
        ];
        assert_eq!(basis, expected_basis);
        Ok(())
    }

    #[rstest]
    #[case("A", vec![1])]
    #[case("B", vec![1])]
    #[case("AB", vec![1,1])]
    #[case("ABB", vec![2, 1])]
    #[case("AAB", vec![2, 1])]
    #[case("ABBB", vec![3, 1])]
    #[case("ABABB", vec![2, 1, 1, 1])]
    #[case("AABAB", vec![2, 1, 1, 1])]
    fn test_generate_goldberg_partitions(
        #[case] word: &str,
        #[case] expected_partition: Vec<usize>,
    ) -> Result<(), LyndonWordError> {
        let partition = word.parse::<LyndonWord<char>>()?.goldberg();
        assert_eq!(partition, expected_partition);
        Ok(())
    }

    #[rstest]
    #[case("AB", vec!["AB", "B"])]
    #[case("AAB", vec!["AAB", "AB", "B"])]
    #[case("AAAB", vec!["AAAB", "AAB", "AB", "B"])]
    #[case("AABAB", vec!["AABAB"])]
    #[case("ABABB", vec!["ABABB"])]
    fn test_generate_right_factors(
        #[case] word: &str,
        #[case] expected_factors: Vec<&str>,
    ) -> Result<(), LyndonWordError> {
        let lyndon_word = word.parse::<LyndonWord<char>>()?;
        let expected_lyndon_factors = expected_factors
            .into_iter()
            .map(|x| x.parse::<LyndonWord<char>>().unwrap())
            .collect::<Vec<_>>();
        let right_factors = lyndon_word.right_factors();
        assert_eq!(right_factors.len(), expected_lyndon_factors.len());
        for (factor, expected_factor) in right_factors
            .into_iter()
            .zip(expected_lyndon_factors.into_iter())
        {
            assert_eq!(factor, expected_factor);
        }
        Ok(())
    }

    #[test]
    fn test_moebius_function() {
        let mu = moebius_mu(50);
        let expected_mu = vec![
            1, -1, -1, 0, -1, 1, -1, 0, 0, 1, // 10
            -1, 0, -1, 1, 1, 0, -1, 0, -1, 0, // 20
            1, 1, -1, 0, 0, 1, 0, 0, -1, -1, // 30
            -1, 0, 1, 1, 1, 0, -1, 1, 1, 0, // 40
            -1, -1, -1, 0, 0, 1, -1, 0, 0, 0, // 50
        ];
        assert_eq!(mu.len(), expected_mu.len());
        for (term, expected_term) in mu.iter().zip(expected_mu.iter()) {
            assert_eq!(term, expected_term);
        }
    }

    #[test]
    fn test_number_of_lyndon_words_per_degree() {
        let basis = LyndonBasis::<u8>::new(2, Sort::Lexicographical);
        let num_words_per_degree = basis.number_of_words_per_degree(5);
        let expected_num_words_per_degree = [2, 1, 2, 3, 6];
        assert_eq!(
            num_words_per_degree.len(),
            expected_num_words_per_degree.len()
        );
        for (term, expected_term) in num_words_per_degree
            .iter()
            .zip(expected_num_words_per_degree.iter())
        {
            assert_eq!(term, expected_term);
        }
    }
}
