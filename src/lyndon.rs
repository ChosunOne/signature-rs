use std::{
    collections::{HashMap, HashSet},
    fmt::{Debug, Display},
    hash::Hash,
    marker::PhantomData,
    ops::Mul,
    str::FromStr,
};

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
        if N > 26 {
            panic!("Only up to 26 generators are supported for 'char' based generators.");
        }

        let alphabet_letters = [
            'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q',
            'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
        ];

        let mut letters = ['A'; N];
        // Special case N == 2 because I want to
        if N != 2 {
            letters.copy_from_slice(&alphabet_letters[..N]);
        } else {
            letters[0] = 'X';
            letters[1] = 'Y';
        }

        letters
    }
}

#[derive(Default, Debug, Clone, Copy)]
pub struct LyndonBasis<const N: usize, T: Generator> {
    _generator: PhantomData<T>,
}

impl<const N: usize, T: Generator<Letter = T>> LyndonBasis<N, T> {
    pub fn generate_basis(k: usize) -> Vec<LyndonWord<N, T>> {
        let alphabet = T::alphabet::<N>();
        dbg!(&alphabet);
        let letter_index: HashMap<_, _> = alphabet
            .iter()
            .copied()
            .enumerate()
            .map(|(i, v)| (v, i))
            .collect();
        let mut unsorted_basis = Vec::new();
        if k == 0 {
            return unsorted_basis;
        }
        let mut w = vec![];

        loop {
            if w.is_empty() {
                w = vec![alphabet[0]];
            } else {
                *w.last_mut().unwrap() = alphabet[letter_index[w.last().unwrap()] + 1];
            }

            if !w.is_empty() && w.len() <= k && *w.last().unwrap() <= alphabet[N - 1] {
                unsorted_basis
                    .push(LyndonWord::try_from(w.clone()).expect("To make a lyndon word"));

                let m = w.len();
                while w.len() < k as usize {
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
        unsorted_basis.sort_by_key(|word| word.len());
        let mut basis = Vec::with_capacity(unsorted_basis.len());
        let mut sorted_basis_index = HashMap::new();

        for level in 1..=k {
            let mut unsorted_words_by_level = unsorted_basis
                .iter()
                .filter(|word| word.len() == level)
                .collect::<Vec<_>>();

            if level == 1 {
                // Sort lexicographically
                unsorted_words_by_level.sort_by_key(|&w| w);
                for word in unsorted_words_by_level {
                    sorted_basis_index.insert(word.clone(), basis.len());
                    basis.push(word.clone());
                }
                continue;
            }

            // Topological Sort
            unsorted_words_by_level.sort_by_key(|w| {
                let (_, i_dp) = w.factorize();
                sorted_basis_index[&i_dp]
            });
            for word in unsorted_words_by_level {
                sorted_basis_index.insert(word.clone(), basis.len());
                basis.push(word.clone());
            }
        }

        basis
    }
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
        return (
            Self::try_from(v).expect("A factorized lyndon word should produce a lyndon word"),
            Self::try_from(w).expect("A factorized lyndon word should produce a lyndon word"),
        );
    }

    pub fn len(&self) -> usize {
        self.letters.len()
    }

    /// Computes the canonical Goldberg representation of the Lyndon word
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

    pub fn right_factors(&self) -> Vec<Self> {
        let alphabet: [_; N] = T::alphabet();
        let mut factors = vec![];
        factors.push(self.clone());
        if self.len() > 1 {
            let (v, w) = self.factorize();
            if v.letters[0] == alphabet[0] {
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
        Ok(Self::try_from([self.letters, rhs.letters].concat())?)
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
        let distinct_word_letters = HashSet::<T>::from_iter(value.iter().copied());
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
    fn it_makes_a_lyndon_word() -> Result<(), LyndonWordError> {
        let letters = "XXY";
        let word = letters.parse::<LyndonWord<2, char>>()?;
        assert_eq!(&format!("{word}"), "XXY");
        Ok(())
    }

    #[test]
    fn it_factorizes_a_lyndon_word() -> Result<(), LyndonWordError> {
        let letters = "XXXXXY";
        let word = letters.parse::<LyndonWord<2, char>>()?;
        let (v, w) = word.factorize();
        assert_eq!(&format!("{v}"), "X");
        assert_eq!(&format!("{w}"), "XXXXY");
        Ok(())
    }

    #[test]
    fn it_grafts_two_lyndon_words() -> Result<(), LyndonWordError> {
        let letters = "XXY";
        let a = letters.parse::<LyndonWord<2, char>>()?;
        let letters = "XY";
        let b = letters.parse::<LyndonWord<2, char>>()?;
        let c = (a * b)?;
        assert_eq!(&format!("{c}"), "XXYXY");
        Ok(())
    }

    #[test]
    fn it_generates_a_lyndon_basis_1() {
        let basis = LyndonBasis::<5, u8>::generate_basis(5);
        assert_eq!(basis.len(), 829);
        for word in &basis {
            assert!(LyndonWord::<5, u8>::is_lyndon(&word.letters));
        }
    }

    #[test]
    fn it_generates_a_lyndon_basis_2() -> Result<(), LyndonWordError> {
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
    fn it_generates_goldberg_partitions(
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
    #[case("XXYXY", vec![])]
    #[case("XYXYY", vec![])]
    fn it_generates_right_factors(
        #[case] word: &str,
        #[case] expected_factors: Vec<&str>,
    ) -> Result<(), LyndonWordError> {
        let lyndon_word = word.parse::<LyndonWord<2, char>>()?;
        let expected_lyndon_factors = expected_factors
            .into_iter()
            .map(|x| x.parse::<LyndonWord<2, char>>().unwrap())
            .collect::<Vec<_>>();
        let right_factors = lyndon_word.right_factors();
        for (factor, expected_factor) in right_factors
            .into_iter()
            .zip(expected_lyndon_factors.into_iter())
        {
            assert_eq!(factor, expected_factor);
        }
        Ok(())
    }
}
