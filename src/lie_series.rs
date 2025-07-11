use crate::{
    Rational,
    bch::goldberg_coeff,
    lyndon::{Generator, LyndonBasis, LyndonWord},
};

/// A Lie series of words with alphabet size `N` and generators of type `T`
pub struct LieSeries<const N: usize, T: Generator<Letter = T>> {
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
    pub goldberg_coefficients: Vec<Rational>,
}

impl<const N: usize, T: Generator<Letter = T>> LieSeries<N, T> {
    pub fn new(max_degree: usize) -> Self {
        let number_of_words_per_degree =
            LyndonBasis::<N, T>::number_of_words_per_degree(max_degree);
        let words = LyndonBasis::<N, T>::generate_basis(max_degree);
        let mut left_factor = vec![0; words.len()];
        let mut right_factor = vec![0; words.len()];
        let mut word_lengths = vec![0; words.len()];
        let mut index_of_degree = vec![0; max_degree + 1];

        for n in 1..=max_degree {
            index_of_degree[n] = index_of_degree[n - 1] + number_of_words_per_degree[n - 1];
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
        }
        let mut goldberg_coefficients = vec![Rational::default(); words.len()];

        let alphabet = T::alphabet::<N>();
        for (i, word) in words.iter().enumerate() {
            let goldberg_partition = word.goldberg();
            let a_first = word.letters[0] == alphabet[0];
            goldberg_coefficients[i] = goldberg_coeff(goldberg_partition, a_first);
        }

        Self {
            max_degree,
            words,
            left_factor,
            right_factor,
            word_lengths,
            index_of_degree,
            goldberg_coefficients,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_lie_series_construction() {
        let lie_series = LieSeries::<2, char>::new(5);
        assert_eq!(lie_series.max_degree, 5);
        assert_eq!(lie_series.words.len(), 14);
        let expected_left_factor = vec![0, 0, 0, 2, 0, 3, 0, 0, 5, 4, 2, 0, 0, 0];
        for (i, (term, expected_term)) in lie_series
            .left_factor
            .iter()
            .zip(expected_left_factor.iter())
            .enumerate()
        {
            dbg!(i);
            assert_eq!(term, expected_term);
        }
        let expected_right_factor = vec![0, 0, 1, 1, 2, 1, 3, 4, 1, 2, 3, 5, 6, 7];
        for (i, (term, expected_term)) in lie_series
            .right_factor
            .iter()
            .zip(expected_right_factor.iter())
            .enumerate()
        {
            dbg!(i);
            assert_eq!(term, expected_term);
        }
        // ABBBB   -1/720
        // ABABB   -1/120
        // AABBB   1/180
        // AABAB   -1/120
        // AAABB   1/180
        // AAAAB   -1/720
        // ABBB    0/1
        // AABB    1/24
        // AAAB    0/1
        // ABB     1/12
        // AAB     1/12
        // AB      1/2
        // B       1/1
        // A       1/1
        let expected_goldberg_coefficients = vec![
            Rational::new(1, 1),  // A
            Rational::new(1, 1),  // B
            Rational::new(1, 2),  // AB
            Rational::new(1, 12), // ABB
            Rational::new(1, 12),
            Rational::default(),
            Rational::new(1, 24),
            Rational::default(),
            Rational::new(-1, 720),
        ];
        dbg!(&lie_series.goldberg_coefficients);
        todo!();
    }
}
