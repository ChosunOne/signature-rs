use std::{
    collections::{HashMap, HashSet},
    fmt::{Debug, Display},
    hash::Hash,
    marker::PhantomData,
    ops::Mul,
    str::FromStr,
};

#[cfg(feature = "progress")]
use indicatif::{ProgressBar, ProgressStyle};
use thiserror::Error;

pub trait Generator:
    Debug + PartialEq + PartialOrd + Ord + Copy + Clone + Display + Eq + Hash
{
    type Letter: Debug + PartialEq + PartialOrd + Ord + Copy + Clone + Display + Eq + Hash;
    fn alphabet<const N: usize>() -> [Self::Letter; N];
}

impl Generator for u8 {
    type Letter = Self;

    fn alphabet<const N: usize>() -> [Self::Letter; N] {
        let mut letters = [0; N];

        for i in 0..N {
            letters[i] = i as u8;
        }

        letters
    }
}

impl Generator for u16 {
    type Letter = Self;

    fn alphabet<const N: usize>() -> [Self::Letter; N] {
        let mut letters = [0; N];

        for i in 0..N {
            letters[i] = i as u16;
        }

        letters
    }
}

impl Generator for u32 {
    type Letter = Self;

    fn alphabet<const N: usize>() -> [Self::Letter; N] {
        let mut letters = [0; N];

        for i in 0..N {
            letters[i] = i as u32;
        }

        letters
    }
}

impl Generator for u64 {
    type Letter = Self;

    fn alphabet<const N: usize>() -> [Self::Letter; N] {
        let mut letters = [0; N];

        for i in 0..N {
            letters[i] = i as u64;
        }

        letters
    }
}

impl Generator for u128 {
    type Letter = Self;

    fn alphabet<const N: usize>() -> [Self::Letter; N] {
        let mut letters = [0; N];

        for i in 0..N {
            letters[i] = i as u128;
        }

        letters
    }
}

impl Generator for i8 {
    type Letter = Self;

    fn alphabet<const N: usize>() -> [Self::Letter; N] {
        let mut letters = [0; N];

        for i in 0..N {
            letters[i] = i as i8;
        }

        letters
    }
}

impl Generator for i16 {
    type Letter = Self;

    fn alphabet<const N: usize>() -> [Self::Letter; N] {
        let mut letters = [0; N];

        for i in 0..N {
            letters[i] = i as i16;
        }

        letters
    }
}

impl Generator for i32 {
    type Letter = Self;

    fn alphabet<const N: usize>() -> [Self::Letter; N] {
        let mut letters = [0; N];

        for i in 0..N {
            letters[i] = i as i32;
        }

        letters
    }
}

impl Generator for i64 {
    type Letter = Self;

    fn alphabet<const N: usize>() -> [Self::Letter; N] {
        let mut letters = [0; N];

        for i in 0..N {
            letters[i] = i as i64;
        }

        letters
    }
}

impl Generator for i128 {
    type Letter = Self;

    fn alphabet<const N: usize>() -> [Self::Letter; N] {
        let mut letters = [0; N];

        for i in 0..N {
            letters[i] = i as i128;
        }

        letters
    }
}

impl Generator for char {
    type Letter = Self;

    fn alphabet<const N: usize>() -> [Self::Letter; N] {
        assert!(
            (N <= 26),
            "Only up to 26 generators are supported for 'char' based generators."
        );

        let alphabet_letters = [
            'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q',
            'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
        ];

        let mut letters = ['A'; N];
        // Special case N == 2 because I want to
        if N == 2 {
            letters[0] = 'X';
            letters[1] = 'Y';
        } else {
            letters.copy_from_slice(&alphabet_letters[..N]);
        }

        letters
    }
}

#[derive(Debug, Default, Clone, Copy)]
pub struct Lexicographical;
#[derive(Debug, Default, Clone, Copy)]
pub struct Topological;

#[derive(Default, Debug, Clone, Copy)]
pub struct LyndonBasis<const N: usize, T: Generator = u8, U = Topological> {
    _generator: PhantomData<T>,
    _sort: PhantomData<U>,
}

impl<const N: usize, T: Generator<Letter = T>, U> LyndonBasis<N, T, U> {
    fn _generate_basis(max_length: usize) -> Vec<LyndonWord<N, T>> {
        #[cfg(feature = "progress")]
        let style = ProgressStyle::with_template(
            "[{eta_precise}] [{bar:35.green/white}] {pos:>2}/{len:2} {msg}",
        )
        .unwrap()
        .progress_chars("=>-");
        let total_words: usize = Self::number_of_words_per_degree(max_length).iter().sum();

        #[cfg(feature = "progress")]
        let pb = ProgressBar::new(total_words as u64).with_style(style.clone());
        #[cfg(feature = "progress")]
        {
            pb.set_message("Generating Lyndon Basis ");
            pb.tick();
        }

        let alphabet = T::alphabet::<N>();
        let letter_index: HashMap<_, _> = alphabet
            .iter()
            .copied()
            .enumerate()
            .map(|(i, v)| (v, i))
            .collect();
        let mut basis = Vec::with_capacity(total_words);
        if max_length == 0 {
            return basis;
        }
        let mut w = vec![];

        loop {
            if w.is_empty() {
                w = vec![alphabet[0]];
            } else {
                *w.last_mut().unwrap() = alphabet[letter_index[w.last().unwrap()] + 1];
            }

            if !w.is_empty() && w.len() <= max_length && *w.last().unwrap() <= alphabet[N - 1] {
                basis.push(LyndonWord::try_from(w.clone()).expect("To make a lyndon word"));

                #[cfg(feature = "progress")]
                pb.inc(1);

                let m = w.len();
                while w.len() < max_length {
                    w.push(w[w.len() % m]);
                }
            }

            while !w.is_empty() && *w.last().unwrap() >= alphabet[N - 1] {
                w.pop();
            }

            if w.is_empty() {
                break;
            }
        }
        basis.sort_by_key(|word| (word.len(), word.letters.clone()));
        #[cfg(feature = "progress")]
        pb.finish();

        basis
    }

    #[must_use]
    pub fn number_of_words_per_degree(max_degree: usize) -> Vec<usize> {
        let mu = moebius_mu(max_degree);
        let mut words_per_degree = vec![0; max_degree];
        for n in 1..=max_degree {
            let mut d = 1_i64;
            let mut h = 0_i64;
            while d * d < n as i64 {
                let quot = n as i64 / d;
                let rem = n as i64 % d;
                if rem == 0 {
                    h += mu[(d - 1) as usize] * (N as i64).pow(quot as u32)
                        + mu[(quot - 1) as usize] * (N as i64).pow(d as u32);
                }
                d += 1;
            }
            if d * d == n as i64 {
                h += mu[(d - 1) as usize] * (N as i64).pow(d as u32);
            }
            words_per_degree[n - 1] = h as usize / n;
        }

        words_per_degree
    }
}

impl<const N: usize, T: Generator<Letter = T>> LyndonBasis<N, T, Topological> {
    #[must_use]
    pub fn generate_basis(max_length: usize) -> Vec<LyndonWord<N, T>> {
        let basis = Self::_generate_basis(max_length);
        Self::topological_sort(&basis)
    }

    fn topological_sort(basis: &[LyndonWord<N, T>]) -> Vec<LyndonWord<N, T>> {
        #[cfg(feature = "progress")]
        let style = ProgressStyle::with_template(
            "[{eta_precise}] [{bar:35.green/white}] {pos:>2}/{len:2} {msg}",
        )
        .unwrap()
        .progress_chars("=>-");
        let max_length = basis.last().map(|x| x.len()).unwrap_or(0);
        let mut sorted_basis = Vec::with_capacity(basis.len());
        let mut sorted_basis_index = HashMap::new();

        #[cfg(feature = "progress")]
        let pb = ProgressBar::new(max_length as u64).with_style(style.clone());
        #[cfg(feature = "progress")]
        {
            pb.set_message("Performing Topological Sort ");
            pb.tick();
        }

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
                #[cfg(feature = "progress")]
                pb.inc(1);

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
            #[cfg(feature = "progress")]
            pb.inc(1);
        }
        #[cfg(feature = "progress")]
        pb.finish();

        sorted_basis
    }
}

impl<const N: usize, T: Generator<Letter = T>> LyndonBasis<N, T, Lexicographical> {
    #[must_use]
    pub fn generate_basis(max_length: usize) -> Vec<LyndonWord<N, T>> {
        Self::_generate_basis(max_length)
    }
}

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

#[derive(Default, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct LyndonWord<const N: usize, T: Generator<Letter = T>> {
    pub letters: Vec<T>,
}

impl<const N: usize, T: Generator<Letter = T>> Display for LyndonWord<N, T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for letter in &self.letters {
            write!(f, "{letter}")?;
        }
        Ok(())
    }
}

impl<const N: usize, T: Generator<Letter = T>> LyndonWord<N, T> {
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

    #[must_use]
    pub fn factorize(&self) -> (LyndonWord<N, T>, LyndonWord<N, T>) {
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

    #[must_use]
    pub fn len(&self) -> usize {
        self.letters.len()
    }

    /// Computes the canonical Goldberg representation of the Lyndon word
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

    #[must_use]
    pub fn right_factors(&self) -> Vec<Self> {
        let alphabet: [_; N] = T::alphabet();
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

impl<const N: usize, T: Generator<Letter = T>> Mul for LyndonWord<N, T> {
    type Output = Result<Self, LyndonWordError>;

    fn mul(self, rhs: Self) -> Self::Output {
        Self::try_from([self.letters, rhs.letters].concat())
    }
}

#[derive(Error, Debug)]
pub enum LyndonWordError {
    #[error("Word is not a valid Lyndon word")]
    InvalidWord,
    #[error("Word uses letters not in the alphabet")]
    InvalidLetter,
}

impl<const N: usize, T: Generator<Letter = T>> TryFrom<Vec<T>> for LyndonWord<N, T> {
    type Error = LyndonWordError;

    fn try_from(value: Vec<T>) -> Result<Self, Self::Error> {
        let distinct_alphabet_letters = HashSet::<T>::from(T::alphabet::<N>());
        let distinct_word_letters = value.iter().copied().collect::<HashSet<T>>();
        let intersection = distinct_alphabet_letters
            .intersection(&distinct_word_letters)
            .collect::<Vec<_>>();
        if distinct_word_letters.len() > distinct_alphabet_letters.len()
            && distinct_word_letters.len() != intersection.len()
        {
            return Err(LyndonWordError::InvalidLetter);
        }
        if !Self::is_lyndon(&value) {
            return Err(LyndonWordError::InvalidWord);
        }
        Ok(Self { letters: value })
    }
}

impl<const N: usize> FromStr for LyndonWord<N, char> {
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
        let letters = "XXY";
        let word = letters.parse::<LyndonWord<2, char>>()?;
        assert_eq!(&format!("{word}"), "XXY");
        Ok(())
    }

    #[test]
    fn test_factorize_a_lyndon_word() -> Result<(), LyndonWordError> {
        let letters = "XXXXXY";
        let word = letters.parse::<LyndonWord<2, char>>()?;
        let (v, w) = word.factorize();
        assert_eq!(&format!("{v}"), "X");
        assert_eq!(&format!("{w}"), "XXXXY");
        Ok(())
    }

    #[test]
    fn test_graft_two_lyndon_words() -> Result<(), LyndonWordError> {
        let letters = "XXY";
        let a = letters.parse::<LyndonWord<2, char>>()?;
        let letters = "XY";
        let b = letters.parse::<LyndonWord<2, char>>()?;
        let c = (a * b)?;
        assert_eq!(&format!("{c}"), "XXYXY");
        Ok(())
    }

    #[test]
    fn test_generate_a_lyndon_basis_1() {
        let basis = LyndonBasis::<5, u8>::generate_basis(5);
        assert_eq!(basis.len(), 829);
        for word in &basis {
            assert!(LyndonWord::<5, u8>::is_lyndon(&word.letters));
        }
    }

    #[test]
    fn test_generate_a_lyndon_basis_2() -> Result<(), LyndonWordError> {
        let basis = LyndonBasis::<2, u8>::generate_basis(5);
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
    #[case("X", vec![1])]
    #[case("Y", vec![1])]
    #[case("XY", vec![1,1])]
    #[case("XYY", vec![2, 1])]
    #[case("XXY", vec![2, 1])]
    #[case("XYYY", vec![3, 1])]
    #[case("XYXYY", vec![2, 1, 1, 1])]
    #[case("XXYXY", vec![2, 1, 1, 1])]
    fn test_generate_goldberg_partitions(
        #[case] word: &str,
        #[case] expected_partition: Vec<usize>,
    ) -> Result<(), LyndonWordError> {
        let partition = word.parse::<LyndonWord<2, char>>()?.goldberg();
        assert_eq!(partition, expected_partition);
        Ok(())
    }

    #[rstest]
    #[case("XY", vec!["XY", "Y"])]
    #[case("XXY", vec!["XXY", "XY", "Y"])]
    #[case("XXXY", vec!["XXXY", "XXY", "XY", "Y"])]
    #[case("XXYXY", vec!["XXYXY"])]
    #[case("XYXYY", vec!["XYXYY"])]
    fn test_generate_right_factors(
        #[case] word: &str,
        #[case] expected_factors: Vec<&str>,
    ) -> Result<(), LyndonWordError> {
        let lyndon_word = word.parse::<LyndonWord<2, char>>()?;
        let expected_lyndon_factors = expected_factors
            .into_iter()
            .map(|x| x.parse::<LyndonWord<2, char>>().unwrap())
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
        for (i, (term, expected_term)) in mu.iter().zip(expected_mu.iter()).enumerate() {
            dbg!(i);
            assert_eq!(term, expected_term);
        }
    }

    #[test]
    fn test_number_of_lyndon_words_per_degree() {
        let num_words_per_degree = LyndonBasis::<2, char>::number_of_words_per_degree(5);
        let expected_num_words_per_degree = [2, 1, 2, 3, 6];
        assert_eq!(
            num_words_per_degree.len(),
            expected_num_words_per_degree.len()
        );
        for (i, (term, expected_term)) in num_words_per_degree
            .iter()
            .zip(expected_num_words_per_degree.iter())
            .enumerate()
        {
            dbg!(i);
            assert_eq!(term, expected_term);
        }
    }
}
