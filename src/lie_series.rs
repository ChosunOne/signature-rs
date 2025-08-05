use std::fmt::Debug;
use std::ops::{Index, IndexMut};

use num_rational::Ratio;

use crate::{
    Commutator, Int,
    lyndon::{Generator, LyndonWord},
};

pub struct LieSeries<const N: usize, T: Generator<Letter = T>, U: Int + Send + Sync> {
    /// The Lyndon basis for the series
    basis: Vec<LyndonWord<N, T>>,
    /// The coefficients for each of the terms in the series
    coefficients: Vec<Ratio<U>>,
    /// The left factors of each word in the basis
    left_factors: Vec<usize>,
    /// The right factors of each word in the basis
    right_factors: Vec<usize>,
    /// The terms of the series
    terms: Vec<Ratio<U>>,
}

impl<const N: usize, T: Generator<Letter = T>, U: Int + Send + Sync> Index<usize>
    for LieSeries<N, T, U>
{
    type Output = Ratio<U>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.terms[index]
    }
}

impl<const N: usize, T: Generator<Letter = T>, U: Int + Send + Sync> IndexMut<usize>
    for LieSeries<N, T, U>
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.terms[index]
    }
}

impl<const N: usize, T: Generator<Letter = T>, U: Int + Send + Sync> LieSeries<N, T, U> {
    pub fn new(
        basis: Vec<LyndonWord<N, T>>,
        coefficients: Vec<Ratio<U>>,
        left_factors: Vec<usize>,
        right_factors: Vec<usize>,
        terms: Vec<Ratio<U>>,
    ) -> Self {
        Self {
            basis,
            coefficients,
            left_factors,
            right_factors,
            terms,
        }
    }
}

impl<const N: usize, T: Generator<Letter = T>, U: Int + Send + Sync> Commutator<&Self>
    for LieSeries<N, T, U>
{
    type Output = Self;

    // For a series this is the following:
    // A: e_1 + 2e_2 + 3[e_1, e_2]
    // B: 4e_1 + 5e_2 + 6[e_1, e_2]
    //
    // [A, B] = [e_1, 4e_1] + [e_1, 5e_2] + [e_1, 6[e_1, e_2]]
    //          + [2e_2, 4e_1] + [2e_2, 5e_2] + [2e_2, 6[e_1, e_2]]
    //          + [3[e_1, e_2], 4e_1] + [3[e_1, e_2], 5e_2] + [3[e_1, e_2], 6[e_1, e_2]]
    //
    //        = [e_1, 5e_2] + [e_1, 6[e_1, e_2]]
    //          + [2e_2, 4e_1] + [2e_2, 6[e_1, e_2]]
    //          + [3[e_1, e_2], 4e_1] + [3[e_1, e_2], 5e_2]
    //
    //        = e_1*5e_2 - 5e_2*e_1 + e_1*6(e_1*e_2 - e_2*e_1) - 6(e_1*e_2 - e_2*e_1)*e_1
    //          + 2e_2*4e_1 - 4e_1*2e_2 + 2e_2 * 6(e_1*e_2 - e_2*e_1) - 6(e_1*e_2 - e_2*e_1)*2e_2
    //          + 3(e_1*e_2 - e_2*e_1)*4e_1 - 4e_1*3(e_1*e_2 - e_2*e_1)
    //          + 3(e_1*e_2 - e_2*e_1)*5e_2 - 5e_2*3(e_1*e_2 - e_2*e_1)
    //
    //        = 5e_1e_2 - 5e_2e_1 + 6e_1e_1e_2 - 6e_1e_2e_1 - 6e_1e_2e_1 + 6e_2e_1e_1
    //          + 8e_2e_1 - 8e_1e_2 + 12e_2e_1e_2 - 12e_2e_2e_1 - 12e_1e_2e_2 + 12e_2e_1e_2
    //          + 12e_1e_2e_1 - 12e_2e_1e_1 - 12e_1e_1e_2 + 12e_1e_2e_1
    //          + 15e_1e_2e_2 - 15e_2e_1e_2 - 15e_2e_1e_2 + 15e_2e_2e_1
    //
    //      Notational shortcut e_1e_2 = e_12
    //        = 5e_12 - 5e_21 + 6e_112 - 6e_121 - 6e_121 + 6e_211
    //          + 8e_21 - 8e_12 + 12e_212 - 12e_221 - 12e_122 + 12e_212
    //          + 12e_121 - 12e_211 - 12e_112 + 12e_121
    //          + 15e_122 - 15e_212 - 15e_212 + 15e_221
    //
    //        = 5e_12 - 8e_12
    //          + 8e_21 - 5e_21
    //          + 6e_112 - 12e_112
    //          + 12e_121 + 12e_121 - 6e_121 - 6e_121
    //          + 15e_122 - 12e_122
    //          + 6e_211 - 12e_211
    //          + 12e_212 + 12e_212 - 15e_212 - 15e_212
    //          + 15e_221 - 12e_221
    //
    //        = -3e_12 + 3e_21 - 6e_112 + 12e_121 + 3e_122 - 6e_211 - 6e_212 + 3e_221
    //
    //  Essentially, for two series A and B with bases:
    //  e_1 + e_2 + [e_1, e_2] + [e_1, [e_1, e_2]] + ...
    //
    //  we need to calculate all the possible combination terms:
    //  [e_1, e_1] + [e_1, e_2] + [e_1, [e_1, [e_1, e_2]]] + ...
    //  + [e_2, e_1] + [e_2, e_2] + [e_2, [e_1, e_2]] + ...
    //
    //  once we have that, we then need to store the corresponding lyndon basis element for each
    //  possible commutation.

    fn commutator(&self, other: &Self) -> Self::Output {
        let mut terms = vec![U::default(); self.terms.len()];
        for i in 0..self.terms.len() {
            for j in 0..other.terms.len() {}
        }
        todo!()
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum CommutatorTerm<T: Debug + Int, U: Clone + Debug + PartialEq + Eq> {
    Atom(U),
    Expression(CommutatorExpression<T, U>),
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum CommutatorExpressionError {
    TooShort,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct CommutatorExpression<T: Debug + Int + PartialEq + Eq, U: Clone + Debug + PartialEq + Eq> {
    coefficient: T,
    left: Box<CommutatorTerm<T, U>>,
    right: Box<CommutatorTerm<T, U>>,
}

impl<T: Debug + Int, U: Clone + Debug + PartialEq + Eq> CommutatorExpression<T, U> {
    /// Implements `[A, B] = -[B, A]`
    pub fn anti_symm(&self) -> Self {
        Self {
            coefficient: -self.coefficient.clone(),
            left: self.right.clone(),
            right: self.left.clone(),
        }
    }
}

impl<const N: usize, T: Int, U: Generator<Letter = U>> TryFrom<&LyndonWord<N, U>>
    for CommutatorExpression<T, U>
{
    type Error = CommutatorExpressionError;

    fn try_from(value: &LyndonWord<N, U>) -> Result<Self, Self::Error> {
        if value.len() <= 1 {
            return Err(CommutatorExpressionError::TooShort);
        }

        let (left, right) = value.factorize();
        let left = if left.len() == 1 {
            Box::new(CommutatorTerm::Atom(left.letters[0]))
        } else {
            Box::new(CommutatorTerm::Expression(Self::try_from(&left)?))
        };
        let right = if right.len() == 1 {
            Box::new(CommutatorTerm::Atom(right.letters[0]))
        } else {
            Box::new(CommutatorTerm::Expression(Self::try_from(&right)?))
        };

        Ok(Self {
            coefficient: T::from(1),
            left,
            right,
        })
    }
}

#[cfg(test)]
mod test {
    use crate::lyndon::LyndonWordError;

    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case("XY", CommutatorExpression::<i128, char> { coefficient: 1, left: Box::new(CommutatorTerm::<i128, char>::Atom('X')), right: Box::new(CommutatorTerm::<i128, char>::Atom('Y'))})]
    fn test_commutator_expression(
        #[case] word: &str,
        #[case] expected_expression: CommutatorExpression<i128, char>,
    ) -> Result<(), LyndonWordError> {
        let word = word.parse::<LyndonWord<2, char>>()?;
        let expression = CommutatorExpression::<i128, char>::try_from(&word)
            .expect("Failed to create expression");
        assert_eq!(expression, expected_expression);

        Ok(())
    }
}
