use std::{
    collections::HashMap,
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
    fn alphabet(size: usize) -> Vec<Self>;
}

impl Generator for u8 {
    fn alphabet(size: usize) -> Vec<Self> {
        let mut letters = vec![0; size];

        for i in 0..size {
            letters[i] = i as u8;
        }

        letters
    }
}

impl Generator for u16 {
    fn alphabet(size: usize) -> Vec<Self> {
        let mut letters = vec![0; size];

        for i in 0..size {
            letters[i] = i as u16;
        }

        letters
    }
}

impl Generator for u32 {
    fn alphabet(size: usize) -> Vec<Self> {
        let mut letters = vec![0; size];

        for i in 0..size {
            letters[i] = i as u32;
        }

        letters
    }
}

impl Generator for u64 {
    fn alphabet(size: usize) -> Vec<Self> {
        let mut letters = vec![0; size];

        for i in 0..size {
            letters[i] = i as u64;
        }

        letters
    }
}

impl Generator for u128 {
    fn alphabet(size: usize) -> Vec<Self> {
        let mut letters = vec![0; size];

        for i in 0..size {
            letters[i] = i as u128;
        }

        letters
    }
}

impl Generator for i8 {
    fn alphabet(size: usize) -> Vec<Self> {
        let mut letters = vec![0; size];

        for i in 0..size {
            letters[i] = i as i8;
        }

        letters
    }
}

impl Generator for i16 {
    fn alphabet(size: usize) -> Vec<Self> {
        let mut letters = vec![0; size];

        for i in 0..size {
            letters[i] = i as i16;
        }

        letters
    }
}

impl Generator for i32 {
    fn alphabet(size: usize) -> Vec<Self> {
        let mut letters = vec![0; size];

        for i in 0..size {
            letters[i] = i as i32;
        }

        letters
    }
}

impl Generator for i64 {
    fn alphabet(size: usize) -> Vec<Self> {
        let mut letters = vec![0; size];

        for i in 0..size {
            letters[i] = i as i64;
        }

        letters
    }
}

impl Generator for i128 {
    fn alphabet(size: usize) -> Vec<Self> {
        let mut letters = vec![0; size];

        for i in 0..size {
            letters[i] = i as i128;
        }

        letters
    }
}

impl Generator for char {
    fn alphabet(size: usize) -> Vec<Self> {
        assert!(
            (size <= 26),
            "Only up to 26 generators are supported for 'char' based generators."
        );

        let alphabet_letters = [
            'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q',
            'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
        ];

        let mut letters = vec!['A'; size];
        letters.copy_from_slice(&alphabet_letters[..size]);

        letters
    }
}

#[derive(Copy, Clone, Default, Debug)]
pub enum Sort {
    #[default]
    Lexicographical,
    Topological,
}

#[derive(Default, Debug, Clone, Copy)]
pub struct LyndonBasis<T: Generator = u8> {
    pub alphabet_size: usize,
    sort: Sort,
    _generator: PhantomData<T>,
}

impl<T: Generator> LyndonBasis<T> {
    pub fn new(alphabet_size: usize, sort: Sort) -> Self {
        Self {
            alphabet_size,
            sort,
            _generator: PhantomData,
        }
    }

    fn _generate_basis(&self, max_length: usize) -> Vec<LyndonWord<T>> {
        #[cfg(feature = "progress")]
        let style = ProgressStyle::with_template(
            "[{eta_precise}] [{bar:35.green/white}] {pos:>2}/{len:2} {msg}",
        )
        .unwrap()
        .progress_chars("=>-");
        let total_words: usize = self.number_of_words_per_degree(max_length).iter().sum();

        #[cfg(feature = "progress")]
        let pb = ProgressBar::new(total_words as u64).with_style(style.clone());
        #[cfg(feature = "progress")]
        {
            pb.set_message("Generating Lyndon Basis ");
            pb.tick();
        }

        let alphabet = T::alphabet(self.alphabet_size);
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

            if !w.is_empty()
                && w.len() <= max_length
                && *w.last().unwrap() <= alphabet[self.alphabet_size - 1]
            {
                basis.push(LyndonWord::try_from(w.clone()).expect("To make a lyndon word"));

                #[cfg(feature = "progress")]
                pb.inc(1);

                let m = w.len();
                while w.len() < max_length {
                    w.push(w[w.len() % m]);
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
        #[cfg(feature = "progress")]
        pb.finish();

        basis
    }

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

    #[must_use]
    pub fn generate_basis(&self, max_length: usize) -> Vec<LyndonWord<T>> {
        let basis = self._generate_basis(max_length);
        match self.sort {
            Sort::Lexicographical => basis,
            Sort::Topological => Self::topological_sort(&basis),
        }
    }

    fn topological_sort(basis: &[LyndonWord<T>]) -> Vec<LyndonWord<T>> {
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
pub struct LyndonWord<T: Generator> {
    pub letters: Vec<T>,
}

impl<T: Generator> Display for LyndonWord<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for letter in &self.letters {
            write!(f, "{letter}")?;
        }
        Ok(())
    }
}

impl<T: Generator> LyndonWord<T> {
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

impl<T: Generator> Mul for LyndonWord<T> {
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

impl<T: Generator> TryFrom<Vec<T>> for LyndonWord<T> {
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
        for (i, (term, expected_term)) in mu.iter().zip(expected_mu.iter()).enumerate() {
            dbg!(i);
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
