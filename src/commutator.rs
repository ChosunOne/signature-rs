use crate::{
    Int,
    lyndon::{Generator, LyndonWord},
};
use std::{fmt::Debug, ops::Mul};

pub trait Commutator<Rhs = Self> {
    type Output;
    /// The commutator operation is represented with `[A, B]` and commonly represents `AB - BA`
    fn commutator(&self, other: Rhs) -> Self::Output;
}

impl<T: Int> Commutator<&Self> for T {
    type Output = T;

    fn commutator(&self, other: &Self) -> Self::Output {
        self.clone() * other.clone() - other.clone() * self.clone()
    }
}

/// Shorthand for applying the commutator operation `[A, B]`.
#[macro_export]
macro_rules! comm {
    ($a:expr, $b:expr) => {
        $a.commutator(&$b)
    };
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum CommutatorTerm<T: Int, U: Clone + Debug + PartialEq + Eq + PartialOrd + Ord> {
    Atom(U),
    Expression(CommutatorExpression<T, U>),
}

impl<T: Debug + Int, U: Clone + Debug + PartialEq + Eq + PartialOrd + Ord> CommutatorTerm<T, U> {
    pub fn lyndon_sort(&mut self) {
        match self {
            // Do nothing, already sorted
            CommutatorTerm::Atom(_) => {}
            CommutatorTerm::Expression(e) => {
                e.left.lyndon_sort();
                e.right.lyndon_sort();
                match e.left.partial_cmp(&e.right) {
                    Some(o) => match o {
                        std::cmp::Ordering::Equal => e.coefficient = T::from(0),
                        std::cmp::Ordering::Greater => {
                            e.coefficient = -e.coefficient.clone();
                            std::mem::swap(&mut e.left, &mut e.right);
                        }
                        _ => {}
                    },
                    None => {}
                }
                // Propagate up coefficients
                if let CommutatorTerm::Expression(e2) = &mut *e.left {
                    e.coefficient *= e2.coefficient.clone();
                    if e2.coefficient == T::from(-1) {
                        e2.coefficient = -e2.coefficient.clone();
                    }
                }
                if let CommutatorTerm::Expression(e2) = &mut *e.right {
                    e.coefficient *= e2.coefficient.clone();
                    if e2.coefficient == T::from(-1) {
                        e2.coefficient = -e2.coefficient.clone();
                    }
                }
            }
        }
    }
}

impl<T: Debug + Int, U: Clone + Debug + PartialEq + Eq + PartialOrd + Ord> Commutator<&Self>
    for CommutatorTerm<T, U>
{
    type Output = Self;

    fn commutator(&self, other: &Self) -> Self::Output {
        match (self, other) {
            (a @ CommutatorTerm::Atom(e1), b @ CommutatorTerm::Atom(e2)) => {
                let coefficient = if e1 == e2 { T::from(0) } else { T::from(1) };
                CommutatorTerm::Expression(CommutatorExpression {
                    coefficient,
                    left: Box::new(a.clone()),
                    right: Box::new(b.clone()),
                })
            }
            (a @ CommutatorTerm::Atom(_), CommutatorTerm::Expression(e)) => {
                let coefficient = e.coefficient.clone();
                let mut right = e.clone();
                right.coefficient = T::from(1);
                let right = Box::new(CommutatorTerm::Expression(right));

                CommutatorTerm::Expression(CommutatorExpression {
                    coefficient,
                    left: Box::new(a.clone()),
                    right,
                })
            }
            (CommutatorTerm::Expression(e), b @ CommutatorTerm::Atom(_)) => {
                let coefficient = e.coefficient.clone();
                let mut left = e.clone();
                left.coefficient = T::from(1);
                let left = Box::new(CommutatorTerm::Expression(left));
                CommutatorTerm::Expression(CommutatorExpression {
                    coefficient,
                    left,
                    right: Box::new(b.clone()),
                })
            }
            (CommutatorTerm::Expression(e1), CommutatorTerm::Expression(e2)) => {
                let coefficient = if e1 != e2 {
                    e1.coefficient.clone() * e2.coefficient.clone()
                } else {
                    T::from(0)
                };
                let mut left = e1.clone();
                left.coefficient = T::from(1);
                let left = Box::new(CommutatorTerm::Expression(left));
                let mut right = e2.clone();
                right.coefficient = T::from(1);
                let right = Box::new(CommutatorTerm::Expression(right));
                CommutatorTerm::Expression(CommutatorExpression {
                    coefficient,
                    left,
                    right,
                })
            }
        }
    }
}

impl<T: Debug + Int, U: Clone + Debug + PartialEq + Eq + PartialOrd + Ord> PartialOrd
    for CommutatorTerm<T, U>
{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        match (self, other) {
            (CommutatorTerm::Atom(a), CommutatorTerm::Atom(b)) => a.partial_cmp(b),
            (CommutatorTerm::Atom(_), CommutatorTerm::Expression(b)) => {
                match self.partial_cmp(&(*b.left)) {
                    Some(o) => match o {
                        std::cmp::Ordering::Equal => Some(std::cmp::Ordering::Less),
                        _ => Some(o),
                    },
                    None => None,
                }
            }
            (CommutatorTerm::Expression(e), CommutatorTerm::Atom(_)) => {
                match (*e.left).partial_cmp(other) {
                    Some(o) => match o {
                        std::cmp::Ordering::Equal => Some(std::cmp::Ordering::Greater),
                        _ => Some(o),
                    },
                    None => None,
                }
            }
            (CommutatorTerm::Expression(s_e), CommutatorTerm::Expression(o_e)) => {
                match s_e.left.partial_cmp(&o_e.left) {
                    Some(o) => match o {
                        std::cmp::Ordering::Equal => s_e.right.partial_cmp(&o_e.right),
                        _ => Some(o),
                    },
                    None => None,
                }
            }
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CommutatorExpressionError {
    TooShort,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct CommutatorExpression<
    T: Debug + Int + PartialEq + Eq,
    U: Clone + Debug + PartialEq + Eq + PartialOrd + Ord,
> {
    pub coefficient: T,
    left: Box<CommutatorTerm<T, U>>,
    right: Box<CommutatorTerm<T, U>>,
}

impl<T: Debug + Int, U: Clone + Debug + PartialEq + Eq + Ord + PartialOrd>
    CommutatorExpression<T, U>
{
    /// Implements `[A, B] = -[B, A]`
    pub fn anti_symm(&self) -> Self {
        Self {
            coefficient: -self.coefficient.clone(),
            left: self.right.clone(),
            right: self.left.clone(),
        }
    }
}

impl<T: Debug + Int, U: Clone + Debug + PartialEq + Eq + PartialOrd + Ord> Commutator
    for CommutatorExpression<T, U>
{
    type Output = Self;

    fn commutator(&self, other: Self) -> Self::Output {
        let coefficient = self.coefficient.clone() * other.coefficient.clone();
        let mut left = self.clone();
        left.coefficient = T::from(1);
        let left = Box::new(CommutatorTerm::Expression(left));
        let mut right = other.clone();
        right.coefficient = T::from(1);
        let right = Box::new(CommutatorTerm::Expression(right));
        CommutatorExpression {
            coefficient,
            left,
            right,
        }
    }
}

impl<const N: usize, T: Int, U: Generator<Letter = U> + Debug + Clone> From<&LyndonWord<N, U>>
    for CommutatorTerm<T, U>
{
    fn from(value: &LyndonWord<N, U>) -> Self {
        if value.len() == 1 {
            return CommutatorTerm::Atom(value.letters[0]);
        }
        CommutatorTerm::Expression(value.try_into().expect("Failed to make expression"))
    }
}

impl<const N: usize, T: Int, U: Generator<Letter = U> + Debug + Clone> TryFrom<&LyndonWord<N, U>>
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

/// This can represent terms of the form `αe_1⋅⋅⋅e_n`
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FormalIndeterminate<T: Clone + Debug + Eq + Ord + PartialEq + PartialOrd, U: Int> {
    coefficient: U,
    symbols: Vec<T>,
}

impl<T: Clone + Debug + Eq + Ord + PartialEq + PartialOrd, U: Int> Mul
    for &FormalIndeterminate<T, U>
{
    type Output = FormalIndeterminate<T, U>;

    /// For formal indeterminates, `e_1 ⊗ e_2 = e_1e_2`
    fn mul(self, rhs: Self) -> Self::Output {
        let mut symbols = self.symbols.clone();
        symbols.append(&mut rhs.symbols.clone());
        let coefficient = self.coefficient.clone() * rhs.coefficient.clone();
        FormalIndeterminate {
            coefficient,
            symbols,
        }
    }
}

impl<T: Clone + Debug + Eq + Ord + PartialEq + PartialOrd, U: Int> Mul<U>
    for &FormalIndeterminate<T, U>
{
    type Output = FormalIndeterminate<T, U>;

    fn mul(self, rhs: U) -> Self::Output {
        FormalIndeterminate {
            coefficient: self.coefficient.clone() * rhs,
            symbols: self.symbols.clone()
        }
    }
}


impl<T: Clone + Debug + Eq + Ord + PartialEq + PartialOrd, U: Int> From<&CommutatorTerm<U, T>>
    for Vec<FormalIndeterminate<T, U>>
{
    fn from(value: &CommutatorTerm<U, T>) -> Self {
        match value {
            CommutatorTerm::Atom(a) => vec![FormalIndeterminate::<T, U> {
                coefficient: U::from(1),
                symbols: vec![a.clone()],
            }],
            CommutatorTerm::Expression(e) => {
                let mut left = Self::from(&*e.left);
                match &*e.left {
                    CommutatorTerm::Expression(l) => {
                        for e_i in left.iter_mut() {
                            e_i.coefficient *= l.coefficient.clone();
                        }
                    }
                    _ => {}
                }
                let mut right = Self::from(&*e.right);
                match &*e.right {
                    CommutatorTerm::Expression(r) => {
                        for e_i in right.iter_mut() {
                            e_i.coefficient *= r.coefficient.clone();
                        }
                    }
                    _ => {}
                }

                let mut terms = vec![];
                for l in left.iter() {
                    for r in right.iter() {
                        terms.push(l * r);
                    }
                }

                for r in right.iter() {
                    for l in left.iter() {
                        terms.push(&(r * l) * U::from(-1));
                    }
                }

                terms
            }
        }
    }
}

#[cfg(test)]
mod test {
    use crate::lyndon::LyndonWordError;

    use super::*;
    use num_bigint::BigInt;
    use rstest::rstest;
    use std::cmp::{Ordering, PartialOrd};

    #[test]
    fn test_commutators_int() {
        assert_eq!(comm![1, 2], 0);
    }

    #[test]
    fn test_commutators_bigint() {
        let a = BigInt::from(1);
        let b = BigInt::from(2);

        assert_eq!(comm![a, b], BigInt::from(0));
    }

    #[rstest]
    #[case(
        CommutatorTerm::<i128, char>::Atom('X'),
        CommutatorTerm::<i128, char>::Atom('Y'),
        CommutatorTerm::Expression(CommutatorExpression {
        coefficient: 1,
        left: Box::new(CommutatorTerm::<i128, char>::Atom('X')),
        right: Box::new(CommutatorTerm::<i128, char>::Atom('Y')),
    }))]
    #[case(
        CommutatorTerm::<i128, char>::Atom('Y'),
        CommutatorTerm::<i128, char>::Atom('X'),
        CommutatorTerm::Expression(CommutatorExpression {
            coefficient: 1,
            left: Box::new(CommutatorTerm::<i128, char>::Atom('Y')),
            right: Box::new(CommutatorTerm::<i128, char>::Atom('X')),
        }))]
    #[case(
        CommutatorTerm::<i128, char>::Expression(CommutatorExpression {
            coefficient: 2,
            left: Box::new(CommutatorTerm::<i128, char>::Atom('X')),
            right: Box::new(CommutatorTerm::<i128, char>::Atom('Y')
        )}),
        CommutatorTerm::<i128, char>::Expression(CommutatorExpression {
            coefficient: 3,
            left: Box::new(CommutatorTerm::<i128, char>::Atom('Y')),
            right: Box::new(CommutatorTerm::<i128, char>::Atom('X'))
        }),
        CommutatorTerm::Expression(CommutatorExpression {
            coefficient: 6,
            left: Box::new(CommutatorTerm::<i128, char>::Expression(CommutatorExpression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::<i128, char>::Atom('X')),
                right: Box::new(CommutatorTerm::<i128, char>::Atom('Y')
            )})),
            right: Box::new(CommutatorTerm::<i128, char>::Expression(CommutatorExpression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::<i128, char>::Atom('Y')),
                right: Box::new(CommutatorTerm::<i128, char>::Atom('X'))
            })),
    }))]
    fn test_commutator_terms(
        #[case] a: CommutatorTerm<i128, char>,
        #[case] b: CommutatorTerm<i128, char>,
        #[case] expected_term: CommutatorTerm<i128, char>,
    ) {
        let term = comm![a, b];
        assert_eq!(term, expected_term);
    }

    #[rstest]
    #[case("XY", CommutatorExpression::<i128, char> { coefficient: 1, left: Box::new(CommutatorTerm::<i128, char>::Atom('X')), right: Box::new(CommutatorTerm::<i128, char>::Atom('Y'))})]
    #[case("XXY", CommutatorExpression::<i128, char> {coefficient: 1, left: Box::new(CommutatorTerm::<i128, char>::Atom('X')), right: Box::new(CommutatorTerm::Expression(CommutatorExpression::<i128, char> { coefficient: 1, left: Box::new(CommutatorTerm::<i128, char>::Atom('X')), right: Box::new(CommutatorTerm::Atom('Y')) })) })]
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

    #[rstest]
    #[case("XY", "XYY", Ordering::Less)]
    #[case("XYY", "XY", Ordering::Greater)]
    #[case("XY", "XY", Ordering::Equal)]
    fn test_commutator_expression_ordering(
        #[case] word_1: &str,
        #[case] word_2: &str,
        #[case] expected_ordering: Ordering,
    ) -> Result<(), LyndonWordError> {
        let word_1 = word_1.parse::<LyndonWord<2, char>>()?;
        let word_2 = word_2.parse::<LyndonWord<2, char>>()?;
        let exp_1 = CommutatorTerm::Expression(
            CommutatorExpression::<i128, char>::try_from(&word_1)
                .expect("Failed to create expression"),
        );
        let exp_2 = CommutatorTerm::Expression(
            CommutatorExpression::<i128, char>::try_from(&word_2)
                .expect("Failed to create expression"),
        );
        assert_eq!(exp_1.partial_cmp(&exp_2), Some(expected_ordering));

        Ok(())
    }

    #[rstest]
    #[case(CommutatorTerm::Atom('X'), CommutatorTerm::Atom('Y'), Ordering::Less)]
    #[case(
        CommutatorTerm::Atom('Y'),
        CommutatorTerm::Atom('X'),
        Ordering::Greater
    )]
    #[case(CommutatorTerm::Atom('X'), CommutatorTerm::Atom('X'), Ordering::Equal)]
    #[case(CommutatorTerm::Atom('X'), CommutatorTerm::Atom('X'), Ordering::Equal)]
    #[case(
        CommutatorTerm::Atom('X'),
        CommutatorTerm::Expression(
            CommutatorExpression::<i128, char> {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom('X')),
                right: Box::new(CommutatorTerm::Atom('Y'))
            }),
        Ordering::Less)]
    #[case(
        CommutatorTerm::Atom('Y'),
        CommutatorTerm::Expression(
            CommutatorExpression::<i128, char> {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom('X')),
                right: Box::new(CommutatorTerm::Atom('Y'))
            }),
        Ordering::Greater)]
    #[case(
        CommutatorTerm::Expression(
            CommutatorExpression::<i128, char> {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom('X')),
                right: Box::new(CommutatorTerm::Atom('Y'))
            }),
        CommutatorTerm::Atom('X'),
        Ordering::Greater)]
    #[case(
        CommutatorTerm::Expression(
            CommutatorExpression::<i128, char> {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom('X')),
                right: Box::new(CommutatorTerm::Atom('Y'))
            }),
        CommutatorTerm::Atom('Y'),
        Ordering::Less)]
    #[case(
        CommutatorTerm::Expression(
            CommutatorExpression::<i128, char> {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom('X')),
                right: Box::new(CommutatorTerm::Atom('Y'))
            }),
        CommutatorTerm::Atom('X'),
        Ordering::Greater)]
    fn test_commutator_term_ordering(
        #[case] term_1: CommutatorTerm<i128, char>,
        #[case] term_2: CommutatorTerm<i128, char>,
        #[case] expected_ordering: Ordering,
    ) {
        assert_eq!(term_1.partial_cmp(&term_2), Some(expected_ordering));
    }

    #[rstest]
    #[case(
        CommutatorTerm::Expression(
            CommutatorExpression::<i128, char> {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom('Y')),
                right: Box::new(CommutatorTerm::Atom('X'))
            }),
        CommutatorTerm::Expression(
            CommutatorExpression::<i128, char> {
                coefficient: -1,
                left: Box::new(CommutatorTerm::Atom('X')),
                right: Box::new(CommutatorTerm::Atom('Y'))
            }),
    )]
    #[case(
        CommutatorTerm::Expression(
            CommutatorExpression::<i128, char> {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom('X')),
                right: Box::new(CommutatorTerm::Atom('Y'))
            }),
        CommutatorTerm::Expression(
            CommutatorExpression::<i128, char> {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom('X')),
                right: Box::new(CommutatorTerm::Atom('Y'))
            }),
    )]
    #[case(
        CommutatorTerm::Expression(
            CommutatorExpression::<i128, char> {
                coefficient: 1,
                left: Box::new(
            CommutatorTerm::Expression(CommutatorExpression::<i128, char> {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom('Y')),
                right: Box::new(CommutatorTerm::Atom('X'))
            })),
                right: Box::new(CommutatorTerm::Atom('Y'))
            }),
        CommutatorTerm::Expression(
            CommutatorExpression::<i128, char> {
                coefficient: -1,
                left: Box::new(
            CommutatorTerm::Expression(CommutatorExpression::<i128, char> {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom('X')),
                right: Box::new(CommutatorTerm::Atom('Y'))
            })),
                right: Box::new(CommutatorTerm::Atom('Y'))
            }),
    )]
    #[case(
        CommutatorTerm::Expression(
            CommutatorExpression::<i128, char> {
                coefficient: 1,
                left: Box::new(
            CommutatorTerm::Expression(CommutatorExpression::<i128, char> {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom('Y')),
                right: Box::new(CommutatorTerm::Atom('X'))
            })),
                right: Box::new(CommutatorTerm::Atom('X'))
            }),
        CommutatorTerm::Expression(
            CommutatorExpression::<i128, char> {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom('X')),
                right: Box::new(
            CommutatorTerm::Expression(CommutatorExpression::<i128, char> {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom('X')),
                right: Box::new(CommutatorTerm::Atom('Y'))
            })),
            }),
    )]
    #[case(
        CommutatorTerm::Expression(
            CommutatorExpression::<i128, char> {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom('X')),
                right: Box::new(CommutatorTerm::Atom('X'))
            }),
        CommutatorTerm::Expression(
            CommutatorExpression::<i128, char> {
                coefficient: 0,
                left: Box::new(CommutatorTerm::Atom('X')),
                right: Box::new(CommutatorTerm::Atom('X'))
            }),
    )]
    #[case(
        CommutatorTerm::Expression(
            CommutatorExpression::<i128, char> {
                coefficient: 1,
                left: Box::new(
            CommutatorTerm::Expression(CommutatorExpression::<i128, char> {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom('Y')),
                right: Box::new(CommutatorTerm::Atom('X'))
            })),
                right: Box::new(
            CommutatorTerm::Expression(CommutatorExpression::<i128, char> {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom('X')),
                right: Box::new(CommutatorTerm::Atom('Y'))
            }))
            }),
        CommutatorTerm::Expression(
            CommutatorExpression::<i128, char> {
                coefficient: 0,
                left: Box::new(
            CommutatorTerm::Expression(CommutatorExpression::<i128, char> {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom('X')),
                right: Box::new(CommutatorTerm::Atom('Y'))
            })),
                right: Box::new(
            CommutatorTerm::Expression(CommutatorExpression::<i128, char> {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom('X')),
                right: Box::new(CommutatorTerm::Atom('Y'))
            })),
            }),
    )]
    fn test_commutator_term_lyndon_sorting(
        #[case] mut term: CommutatorTerm<i128, char>,
        #[case] expected_term: CommutatorTerm<i128, char>,
    ) {
        term.lyndon_sort();
        assert_eq!(term, expected_term);
    }

    #[rstest]
    #[case(
        CommutatorTerm::Atom('X'),
        CommutatorTerm::Atom('Y'),
        vec![
            FormalIndeterminate::<char, i128> { 
                coefficient: 1, 
                symbols: vec!['X', 'Y']
            },
            FormalIndeterminate::<char, i128> {
                coefficient: -1,
                symbols: vec!['Y', 'X']
            }
        ])]
    #[case(
        comm![CommutatorTerm::Atom('A'), CommutatorTerm::Atom('B')],
        comm![CommutatorTerm::Atom('C'), CommutatorTerm::Atom('D')],
        vec![
            FormalIndeterminate::<char, i128> {
                coefficient: 1,
                symbols: vec!['A', 'B', 'C', 'D']
            },
            FormalIndeterminate::<char, i128> {
                coefficient: -1,
                symbols: vec!['A', 'B', 'D', 'C']
            },
            FormalIndeterminate::<char, i128> {
                coefficient: -1,
                symbols: vec!['B', 'A', 'C', 'D']
            },
            FormalIndeterminate::<char, i128> {
                coefficient: 1,
                symbols: vec!['B', 'A', 'D', 'C']
            },
            FormalIndeterminate::<char, i128> {
                coefficient: -1,
                symbols: vec!['C', 'D', 'A', 'B']
            },
            FormalIndeterminate::<char, i128> {
                coefficient: 1,
                symbols: vec!['C', 'D', 'B', 'A']
            },
            FormalIndeterminate::<char, i128> {
                coefficient: 1,
                symbols: vec!['D', 'C', 'A', 'B']
            },
            FormalIndeterminate::<char, i128> {
                coefficient: -1,
                symbols: vec!['D', 'C', 'B', 'A']
            }
        ]
    )]
    fn test_formal_indeterminate_expansion(
        #[case] a: CommutatorTerm<i128, char>,
        #[case] b: CommutatorTerm<i128, char>,
        #[case] expected_indeterminates: Vec<FormalIndeterminate<char, i128>>,
    ) {
        let term = comm![a, b];
        let indeterminates = Vec::<FormalIndeterminate<char, i128>>::from(&term);
        assert_eq!(indeterminates, expected_indeterminates);
    }
}
