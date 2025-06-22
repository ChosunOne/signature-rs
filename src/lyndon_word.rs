use std::{collections::HashMap, fmt::Display, hash::Hash, ops::Mul, str::FromStr};

use thiserror::Error;

#[derive(Default, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct LyndonWord<T: Ord + Copy + Clone + Display + Eq + Hash> {
    pub letters: Vec<T>,
}

impl<T: Ord + Copy + Clone + Display + Eq + Hash> Display for LyndonWord<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for letter in &self.letters {
            write!(f, "{letter}")?;
        }
        Ok(())
    }
}

impl<T: Ord + Copy + Clone + Display + Eq + Hash> LyndonWord<T> {
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
        return (
            Self::try_from(v).expect("A factorized lyndon word should produce a lyndon word"),
            Self::try_from(w).expect("A factorized lyndon word should produce a lyndon word"),
        );
    }

    pub fn len(&self) -> usize {
        self.letters.len()
    }
}

impl<T: Ord + Copy + Clone + Display + Hash> Mul for LyndonWord<T> {
    type Output = Result<Self, LyndonWordError>;

    fn mul(self, rhs: Self) -> Self::Output {
        Ok(Self::try_from([self.letters, rhs.letters].concat())?)
    }
}

impl LyndonWord<u8> {
    pub fn generate_basis(n: u8, k: usize) -> Vec<LyndonWord<u8>> {
        let mut unsorted_basis = Vec::new();
        if k == 0 {
            return unsorted_basis;
        }
        let mut w = vec![];

        loop {
            if w.is_empty() {
                w = vec![0_u8];
            } else {
                *w.last_mut().unwrap() += 1;
            }

            if !w.is_empty() && w.len() <= k && *w.last().unwrap() < n {
                unsorted_basis
                    .push(LyndonWord::try_from(w.clone()).expect("To make a lyndon word"));

                let m = w.len();
                while w.len() < k as usize {
                    w.push(w[w.len() % m]);
                }
            }

            while !w.is_empty() && *w.last().unwrap() >= n {
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

#[derive(Error, Debug)]
pub enum LyndonWordError {
    #[error("Word is not a valid Lyndon word")]
    InvalidWord,
}

impl<T: Ord + Copy + Clone + Display + Hash> TryFrom<Vec<T>> for LyndonWord<T> {
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

    use super::*;

    #[test]
    fn it_makes_a_lyndon_word() -> Result<(), LyndonWordError> {
        let letters = "xxy";
        let word = letters.parse::<LyndonWord<char>>()?;
        assert_eq!(&format!("{word}"), "xxy");
        Ok(())
    }

    #[test]
    fn it_factorizes_a_lyndon_word() -> Result<(), LyndonWordError> {
        let letters = "xxxxxy";
        let word = letters.parse::<LyndonWord<char>>()?;
        let (v, w) = word.factorize();
        assert_eq!(&format!("{v}"), "x");
        assert_eq!(&format!("{w}"), "xxxxy");
        Ok(())
    }

    #[test]
    fn it_grafts_two_lyndon_words() -> Result<(), LyndonWordError> {
        let letters = "xxy";
        let a = letters.parse::<LyndonWord<char>>()?;
        let letters = "xy";
        let b = letters.parse::<LyndonWord<char>>()?;
        let c = (a * b)?;
        assert_eq!(&format!("{c}"), "xxyxy");
        Ok(())
    }

    #[test]
    fn it_generates_a_lyndon_basis_1() {
        let basis = LyndonWord::<u8>::generate_basis(5, 5);
        assert_eq!(basis.len(), 829);
        for word in &basis {
            assert!(LyndonWord::<u8>::is_lyndon(&word.letters));
        }
    }

    #[test]
    fn it_generates_a_lyndon_basis_2() -> Result<(), LyndonWordError> {
        let basis = LyndonWord::<u8>::generate_basis(2, 5);
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
}
