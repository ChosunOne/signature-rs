use crate::{
    lyndon::{Generator, LyndonWord}, Arith
};
use std::{fmt::{Debug, Display}, hash::Hash, ops::{Mul, Neg}};


pub trait Commutator<Rhs = Self> {
    type Output;
    /// The commutator operation is represented with `[A, B]` and commonly represents `AB - BA`
    fn commutator(&self, other: Rhs) -> Self::Output;
}

impl<T: Arith> Commutator<&Self> for T {
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
pub enum CommutatorTerm<T: Arith, U: Clone + Debug + Hash + PartialEq + PartialOrd> {
    Atom{
        coefficient: T,
        atom: U,
    },
    Expression{
        coefficient: T,
        left: Box<Self>,
        right: Box<Self>,
    },
}

impl<T: Arith + Display, U: Clone + Display + Debug + PartialEq + PartialOrd + Hash> Display for CommutatorTerm<T, U> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Atom{ coefficient, atom} => {
                if coefficient.is_one() {
                    write!(f, "{atom}")
                } else {
                    write!(f, "{coefficient} * {atom}")
                }
            }
            Self::Expression { coefficient, left, right} => {
                if coefficient.is_one() {
                    write!(f, "[{left}, {right}]")
                } else {
                    write!(f, "{coefficient} * [{left}, {right}]")
                }
            }
        }
    }
}

impl<T: Arith, U: Clone + Debug + PartialEq + PartialOrd + Hash> Mul<T> for CommutatorTerm<T, U> {
    type Output = Self;

    fn mul(self, rhs: T) -> Self::Output {
        match self {
            Self::Atom { coefficient, atom } => {
                Self::Atom {
                    coefficient: coefficient * rhs,
                    atom
                }
            } ,
            Self::Expression { coefficient, left, right } => {
                Self::Expression { coefficient: coefficient * rhs, left, right }
            }
        }
    }
}

impl<T: Arith, U: Clone + Debug + PartialEq + PartialOrd + Hash> Neg for CommutatorTerm<T, U> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        match self {
            Self::Atom { coefficient, atom} => {
                Self::Atom {
                    coefficient: coefficient.neg(),
                    atom
                }
            },
            Self::Expression { coefficient, left, right} => {
                Self::Expression {
                    coefficient: coefficient.neg(),
                    left,
                    right
                }
            }
        }
    }
}


impl<T: Arith, U: Clone + Debug + PartialEq + PartialOrd + Hash> CommutatorTerm<T, U> {
    pub fn lyndon_sort(&mut self) {
        match self {
            // Do nothing, already sorted
            CommutatorTerm::Atom { .. } => {}
            CommutatorTerm::Expression{coefficient,  left, right} => {
                left.lyndon_sort();
                right.lyndon_sort();
                if let Some(o) = left.partial_cmp(&right) { match o {
                    std::cmp::Ordering::Equal => *coefficient = T::zero(),
                    std::cmp::Ordering::Greater => {
                        *coefficient = -coefficient.clone();
                        std::mem::swap(left, right);
                    }
                    _ => {}
                } }
                // Propagate up coefficients
                if let CommutatorTerm::Expression{ coefficient: c2, ..} = &mut **left {
                    *coefficient *= c2.clone();
                    if *c2 == -T::one() {
                        *c2 = -c2.clone();
                    }
                }
                if let CommutatorTerm::Expression{ coefficient: c2, ..} = &mut **right {
                    *coefficient *= c2.clone();
                    if *c2 == -T::one() {
                        *c2 = -c2.clone();
                    }
                }
            }
        }
    }

    /// Uses the anti-symmetry property to produce the form `[[A, B], C] = [A, [B, C]] - [B, [A, C]]`
    pub fn jacobi_identity(&self) -> Option<(Self, Self)> {
        match self {
            CommutatorTerm::Atom{ .. } => None,
            CommutatorTerm::Expression { coefficient, left, right} => {
                // Check if expression is already in basis form
                let CommutatorTerm::Expression{ left: l_left, right: l_right, ..} = &**left else { return None};
                let a = l_left.clone();
                let b = l_right.clone();
                let c = right.clone();

                let mut left_term = comm![a, comm![b, c]] * coefficient.clone();
                left_term.lyndon_sort();
                let mut right_term = comm![b, comm![a, c]] * -coefficient.clone();
                right_term.lyndon_sort();

                Some((left_term, right_term))
            },
        }
    }
}

impl<T: Arith, U: Clone + Debug + PartialEq + PartialOrd + Hash> Commutator<&Self>
    for CommutatorTerm<T, U>
{
    type Output = Self;

    fn commutator(&self, other: &Self) -> Self::Output {
        match (self, other) {
            (a @ CommutatorTerm::Atom {coefficient: c1, atom: a1}, 
                b @ CommutatorTerm::Atom {coefficient: c2, atom: a2}) => {
                let coefficient = if  a == b { T::zero() } else { c1.clone() * c2.clone() };
                let left = Box::new(CommutatorTerm::Atom { coefficient: T::one(), atom: a1.clone()});
                let right = Box::new(CommutatorTerm::Atom { coefficient: T::one(), atom: a2.clone()});
                CommutatorTerm::Expression{
                    coefficient,
                    left,
                    right,
                }
            }
            (CommutatorTerm::Atom {coefficient: c1, atom}, 
                CommutatorTerm::Expression {coefficient: c2, left: l1, right}) => {
                let coefficient = c1.clone() * c2.clone();
                let left = Box::new(CommutatorTerm::Atom { coefficient: T::one(), atom: atom.clone()});
                let right = Box::new(CommutatorTerm::Expression { coefficient: T::one(), left: l1.clone(), right: right.clone()});

                CommutatorTerm::Expression{
                    coefficient,
                    left,
                    right,
                }
            }
            (CommutatorTerm::Expression {coefficient: c1, left: l1, right: r1}, 
                CommutatorTerm::Atom{ coefficient: c2, atom }) => {
                let coefficient = c1.clone() * c2.clone();
                let left = Box::new(CommutatorTerm::Expression { coefficient: T::one(), left: l1.clone(), right: r1.clone() });
                let right = Box::new(CommutatorTerm::Atom { coefficient: T::one(), atom: atom.clone()});
                CommutatorTerm::Expression {
                    coefficient,
                    left,
                    right
                }
            }
            (a @ CommutatorTerm::Expression {coefficient: c1, left: l1, right: r1},
                b @ CommutatorTerm::Expression {coefficient: c2, left: l2, right: r2}) => {
                let coefficient = if a == b {
                    T::zero()
                } else {
                    c1.clone() * c2.clone()
                };
                let left = Box::new(CommutatorTerm::Expression { coefficient: T::one(), left: l1.clone(), right: r1.clone()});
                let right = Box::new(CommutatorTerm::Expression { coefficient: T::one(), left: l2.clone(), right: r2.clone()});

                CommutatorTerm::Expression {
                    coefficient,
                    left,
                    right,
                }
            }
        }
    }
}

impl<T: Arith , U: Clone + Debug + PartialEq + PartialOrd + Hash> PartialOrd
    for CommutatorTerm<T, U>
{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        match (self, other) {
            (CommutatorTerm::Atom { atom: a1, ..}, CommutatorTerm::Atom { atom: a2, ..}) => a1.partial_cmp(a2),
            (CommutatorTerm::Atom { .. }, CommutatorTerm::Expression { left, .. }) => {
                match self.partial_cmp(&(*left)) {
                    Some(o) => match o {
                        std::cmp::Ordering::Equal => Some(std::cmp::Ordering::Less),
                        _ => Some(o),
                    },
                    None => None,
                }
            }
            (CommutatorTerm::Expression { left, ..}, CommutatorTerm::Atom { .. }) => {
                match (**left).partial_cmp(other) {
                    Some(o) => match o {
                        std::cmp::Ordering::Equal => Some(std::cmp::Ordering::Greater),
                        _ => Some(o),
                    },
                    None => None,
                }
            }
            (CommutatorTerm::Expression {left: l1, right: r1, ..}, 
                CommutatorTerm::Expression {left: l2, right: r2, ..}) => {
                match l1.partial_cmp(&l2) {
                    Some(o) => match o {
                        std::cmp::Ordering::Equal => r1.partial_cmp(&r2),
                        _ => Some(o),
                    },
                    None => None,
                }
            }
        }
    }
}

impl<T: Arith, U: Generator + Debug + Clone + Hash> From<&LyndonWord<U>>
    for CommutatorTerm<T, U>
{
    fn from(value: &LyndonWord<U>) -> Self {
        if value.len() == 1 {
            return CommutatorTerm::Atom{coefficient: T::one(), atom: value.letters[0]};
        }

        let (left, right) = value.factorize();
        let left = if left.len() == 1 {
            Box::new(CommutatorTerm::Atom {coefficient: T::one(), atom: left.letters[0]})
        } else {
            Box::new(Self::from(&left))
        };
        let right = if right.len() == 1 {
            Box::new(CommutatorTerm::Atom {coefficient: T::one(), atom: right.letters[0]})
        } else {
            Box::new(Self::from(&right))
        };

        Self::Expression { coefficient: T::one(), left, right }
    }
}


/// This can represent terms of the form `αe_1⋅⋅⋅e_n`
#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct FormalIndeterminate<T: Clone + Debug + Eq + Ord + PartialEq + PartialOrd + Hash, U: Arith> {
    coefficient: U,
    symbols: Vec<T>,
}

impl<T: Clone + Debug + Eq + Ord + PartialEq + PartialOrd + Hash, U: Arith> Mul
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

impl<T: Clone + Debug + Eq + Ord + PartialEq + PartialOrd + Hash, U: Arith> Mul<U>
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

impl<T: Clone + Debug + Eq + Ord + PartialEq + PartialOrd + Hash, U: Arith> Neg for FormalIndeterminate<T, U> {
    type Output = Self;

    fn neg(mut self) -> Self::Output {
        self.coefficient = -self.coefficient;
        self
    }
}

impl<T: Clone + Debug + Eq + Ord + PartialEq + PartialOrd + Hash, U: Arith> Neg for &FormalIndeterminate<T, U> {
    type Output = FormalIndeterminate<T, U>;

    fn neg(self) -> Self::Output {
        let mut result = self.clone();
        result.coefficient = -result.coefficient;
        result
    }
}


impl<T: Clone + Debug + Display + Eq + Ord + PartialEq + PartialOrd + Hash, U: Arith + Display> From<&CommutatorTerm<U, T>>
    for Vec<FormalIndeterminate<T, U>>
{
    fn from(value: &CommutatorTerm<U, T>) -> Self {
        match value {
            CommutatorTerm::Atom { coefficient, atom} => vec![FormalIndeterminate::<T, U> {
                coefficient: coefficient.clone(),
                symbols: vec![atom.clone()],
            }],
            CommutatorTerm::Expression { left: l1, right: r1, ..} => {
                let mut left = Self::from(&**l1);
                if let CommutatorTerm::Expression { coefficient: c_l, ..} = &**l1 {
                    for e_i in &mut left {
                        e_i.coefficient *= c_l.clone();
                    }
                }
                let mut right = Self::from(&**r1);
                if let CommutatorTerm::Expression { coefficient: c_r, ..} = &**r1 {
                    for e_i in &mut right {
                        e_i.coefficient *= c_r.clone();
                    }
                }

                let mut terms = vec![];
                for l in &left {
                    for r in &right {
                        terms.push(l * r);
                    }
                }

                for r in &right {
                    for l in &left {
                        terms.push(-&(r * l));
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
        CommutatorTerm::<i128, char>::Atom {atom: 'A', coefficient: 1},
        CommutatorTerm::<i128, char>::Atom {atom: 'B', coefficient: 1},
        CommutatorTerm::Expression {
        coefficient: 1,
        left: Box::new(CommutatorTerm::<i128, char>::Atom { atom: 'A', coefficient: 1}),
        right: Box::new(CommutatorTerm::<i128, char>::Atom { atom: 'B', coefficient: 1}),
    })]
    #[case(
        CommutatorTerm::<i128, char>::Atom {atom: 'B', coefficient: 1},
        CommutatorTerm::<i128, char>::Atom {atom: 'A', coefficient: 1},
        CommutatorTerm::Expression {
            coefficient: 1,
            left: Box::new(CommutatorTerm::<i128, char>::Atom { atom: 'B', coefficient: 1 }),
            right: Box::new(CommutatorTerm::<i128, char>::Atom { atom: 'A', coefficient: 1 }),
        })]
    #[case(
        CommutatorTerm::<i128, char>::Expression {
            coefficient: 2,
            left: Box::new(CommutatorTerm::<i128, char>::Atom { atom: 'A', coefficient: 1 }),
            right: Box::new(CommutatorTerm::<i128, char>::Atom { atom: 'B', coefficient: 1 })
        },
        CommutatorTerm::<i128, char>::Expression {
            coefficient: 3,
            left: Box::new(CommutatorTerm::<i128, char>::Atom { atom: 'B', coefficient: 1}),
            right: Box::new(CommutatorTerm::<i128, char>::Atom { atom: 'A', coefficient: 1})
        },
        CommutatorTerm::Expression {
            coefficient: 6,
            left: Box::new(CommutatorTerm::<i128, char>::Expression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::<i128, char>::Atom { atom: 'A', coefficient: 1}),
                right: Box::new(CommutatorTerm::<i128, char>::Atom { atom: 'B', coefficient: 1})
            }),
            right: Box::new(CommutatorTerm::<i128, char>::Expression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::<i128, char>::Atom { atom: 'B', coefficient: 1}),
                right: Box::new(CommutatorTerm::<i128, char>::Atom { atom: 'A', coefficient: 1})
            }),
    })]
    fn test_commutator_terms(
        #[case] a: CommutatorTerm<i128, char>,
        #[case] b: CommutatorTerm<i128, char>,
        #[case] expected_term: CommutatorTerm<i128, char>,
    ) {
        let term = comm![a, b];
        assert_eq!(term, expected_term);
    }


    #[rstest]
    #[case("AB", "ABB", Ordering::Less)]
    #[case("ABB", "AB", Ordering::Greater)]
    #[case("AB", "AB", Ordering::Equal)]
    fn test_commutator_expression_ordering(
        #[case] word_1: &str,
        #[case] word_2: &str,
        #[case] expected_ordering: Ordering,
    ) -> Result<(), LyndonWordError> {
        let word_1 = word_1.parse::<LyndonWord<char>>()?;
        let word_2 = word_2.parse::<LyndonWord<char>>()?;
        let exp_1 = CommutatorTerm::<i32, char>::from(&word_1);
        let exp_2 = CommutatorTerm::from(&word_2);
        assert_eq!(exp_1.partial_cmp(&exp_2), Some(expected_ordering));

        Ok(())
    }

    #[rstest]
    #[case(CommutatorTerm::Atom { atom: 'A', coefficient: 1 }, CommutatorTerm::Atom { atom: 'B', coefficient: 1 }, Ordering::Less)]
    #[case(
        CommutatorTerm::Atom { atom: 'B', coefficient: 1 },
        CommutatorTerm::Atom { atom: 'A', coefficient: 1 },
        Ordering::Greater
    )]
    #[case(CommutatorTerm::Atom { atom: 'A', coefficient: 1 }, CommutatorTerm::Atom { atom: 'A', coefficient: 1 }, Ordering::Equal)]
    #[case(CommutatorTerm::Atom { atom: 'A', coefficient: 1 }, CommutatorTerm::Atom { atom: 'A', coefficient: 1 }, Ordering::Equal)]
    #[case(
        CommutatorTerm::Atom { atom: 'A', coefficient: 1 },
        CommutatorTerm::Expression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom { atom: 'A', coefficient: 1 }),
                right: Box::new(CommutatorTerm::Atom { atom: 'B', coefficient: 1 })
            },
        Ordering::Less)]
    #[case(
        CommutatorTerm::Atom { atom: 'B', coefficient: 1 },
        CommutatorTerm::Expression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom { atom: 'A', coefficient: 1 }),
                right: Box::new(CommutatorTerm::Atom { atom: 'B', coefficient: 1 })
            },
        Ordering::Greater)]
    #[case(
        CommutatorTerm::Expression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom { atom: 'A', coefficient: 1 }),
                right: Box::new(CommutatorTerm::Atom { atom: 'B', coefficient: 1 })
            },
        CommutatorTerm::Atom { atom: 'A', coefficient: 1 },
        Ordering::Greater)]
    #[case(
        CommutatorTerm::Expression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom { atom: 'A', coefficient: 1 }),
                right: Box::new(CommutatorTerm::Atom { atom: 'B', coefficient: 1 })
            },
        CommutatorTerm::Atom { atom: 'B', coefficient: 1 },
        Ordering::Less)]
    #[case(
        CommutatorTerm::Expression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom { atom: 'A', coefficient: 1 }),
                right: Box::new(CommutatorTerm::Atom { atom: 'B', coefficient: 1 })
            },
        CommutatorTerm::Atom { atom: 'A', coefficient: 1 },
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
        CommutatorTerm::Expression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom { atom: 'B', coefficient: 1 }),
                right: Box::new(CommutatorTerm::Atom { atom: 'A', coefficient: 1 })
            },
        CommutatorTerm::Expression {
                coefficient: -1,
                left: Box::new(CommutatorTerm::Atom { atom: 'A', coefficient: 1 }),
                right: Box::new(CommutatorTerm::Atom { atom: 'B', coefficient: 1 })
            },
    )]
    #[case(
        CommutatorTerm::Expression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom { atom: 'A', coefficient: 1 }),
                right: Box::new(CommutatorTerm::Atom { atom: 'B', coefficient: 1 })
            },
        CommutatorTerm::Expression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom { atom: 'A', coefficient: 1 }),
                right: Box::new(CommutatorTerm::Atom { atom: 'B', coefficient: 1 })
            },
    )]
    #[case(
        CommutatorTerm::Expression {
                coefficient: 1,
                left: Box::new(
            CommutatorTerm::Expression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom { atom: 'B', coefficient: 1 }),
                right: Box::new(CommutatorTerm::Atom { atom: 'A', coefficient: 1 })
            }),
                right: Box::new(CommutatorTerm::Atom { atom: 'B', coefficient: 1 })
            },
        CommutatorTerm::Expression {
                coefficient: -1,
                left: Box::new(
            CommutatorTerm::Expression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom { atom: 'A', coefficient: 1 }),
                right: Box::new(CommutatorTerm::Atom { atom: 'B', coefficient: 1 })
            }),
                right: Box::new(CommutatorTerm::Atom { atom: 'B', coefficient: 1 })
            },
    )]
    #[case(
        CommutatorTerm::Expression { 
                coefficient: 1,
                left: Box::new(
            CommutatorTerm::Expression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom { atom: 'B', coefficient: 1 }),
                right: Box::new(CommutatorTerm::Atom { atom: 'A', coefficient: 1 })
            }),
                right: Box::new(CommutatorTerm::Atom { atom: 'A', coefficient: 1 })
            },
        CommutatorTerm::Expression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom { atom: 'A', coefficient: 1 }),
                right: Box::new(
            CommutatorTerm::Expression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom { atom: 'A', coefficient: 1 }),
                right: Box::new(CommutatorTerm::Atom { atom: 'B', coefficient: 1 })
            }),
            },
    )]
    #[case(
        CommutatorTerm::Expression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom { atom: 'A', coefficient: 1 }),
                right: Box::new(CommutatorTerm::Atom { atom: 'A', coefficient: 1 })
            },
        CommutatorTerm::Expression {
                coefficient: 0,
                left: Box::new(CommutatorTerm::Atom { atom: 'A', coefficient: 1 }),
                right: Box::new(CommutatorTerm::Atom { atom: 'A', coefficient: 1 })
            },
    )]
    #[case(
        CommutatorTerm::Expression {
                coefficient: 1,
                left: Box::new(
            CommutatorTerm::Expression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom { atom: 'B', coefficient: 1 }),
                right: Box::new(CommutatorTerm::Atom { atom: 'A', coefficient: 1 })
            }),
                right: Box::new(
            CommutatorTerm::Expression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom { atom: 'A', coefficient: 1 }),
                right: Box::new(CommutatorTerm::Atom { atom: 'B', coefficient: 1 })
            })
            },
        CommutatorTerm::Expression {
                coefficient: 0,
                left: Box::new(
            CommutatorTerm::Expression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom { atom: 'A', coefficient: 1 }),
                right: Box::new(CommutatorTerm::Atom { atom: 'B', coefficient: 1 })
            }),
                right: Box::new(
            CommutatorTerm::Expression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::Atom { atom: 'A', coefficient: 1 }),
                right: Box::new(CommutatorTerm::Atom { atom: 'B', coefficient: 1 })
            }),
            },
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
        comm![
            comm![
                CommutatorTerm::Atom { atom: 'A', coefficient: 1 },
                CommutatorTerm::Atom { atom: 'B', coefficient: 1 }],
            CommutatorTerm::Atom { atom: 'C', coefficient: 1 }],
        comm![
            CommutatorTerm::Atom { atom: 'A', coefficient: 1 },
            comm![
                CommutatorTerm::Atom { atom: 'B', coefficient: 1 },
                CommutatorTerm::Atom { atom: 'C', coefficient: 1 }]
        ], 
        CommutatorTerm::Expression {
            coefficient: 1,
            left: Box::new(comm![
                CommutatorTerm::Atom { atom: 'A', coefficient: 1 },
                CommutatorTerm::Atom { atom: 'C', coefficient: 1 }]),
            right: Box::new(CommutatorTerm::Atom { atom: 'B', coefficient: 1 }),
    })]
    fn test_commutator_term_jacobi_identity(
        #[case] term: CommutatorTerm<i128, char>,
        #[case] expected_left_term: CommutatorTerm<i128, char>,
        #[case] expected_right_term: CommutatorTerm<i128, char>) {
        let Some((left_term, right_term)) = term.jacobi_identity() else {
            panic!("Failed to create jacobi identity");
        };
        assert_eq!(left_term, expected_left_term, "{left_term} != {expected_left_term}");
        assert_eq!(right_term, expected_right_term, "{right_term} != {expected_right_term}");
    }

    #[rstest]
    #[case(
        CommutatorTerm::Atom { atom: 'A', coefficient: 1 },
        CommutatorTerm::Atom { atom: 'B', coefficient: 1 },
        vec![
            FormalIndeterminate::<char, i128> { 
                coefficient: 1,
                symbols: vec!['A', 'B']
            },
            FormalIndeterminate::<char, i128> {
                coefficient: -1,
                symbols: vec!['B', 'A']
            }
        ])]
    #[case(
        comm![CommutatorTerm::Atom { atom: 'A', coefficient: 1 }, CommutatorTerm::Atom { atom: 'B', coefficient: 1 }],
        comm![CommutatorTerm::Atom { atom: 'C', coefficient: 1 }, CommutatorTerm::Atom { atom: 'D', coefficient: 1 }],
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
