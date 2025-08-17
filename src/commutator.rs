use ordered_float::NotNan;

use crate::{
    generators::Generator,
    lyndon::LyndonWord, Arith
};
use std::{collections::{HashMap, HashSet}, fmt::{Debug, Display}, hash::Hash, ops::{Mul, Neg}};


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
pub enum CommutatorTerm<T: Arith, U: Clone + Debug + Hash + PartialEq + PartialOrd + Ord> {
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

impl<T: Arith> From<char> for CommutatorTerm<T, char> {
    fn from(value: char) -> Self {
        Self::Atom {
            coefficient: T::one(),
            atom: value
        }
    }
}

impl<T: Arith> From<u8> for CommutatorTerm<T, u8> {
    fn from(value: u8) -> Self {
        Self::Atom {
            coefficient: T::one(),
            atom: value
        }
    }
}

impl<T: Arith + Display, U: Clone + Display + Debug + PartialEq + PartialOrd + Ord + Hash> Display for CommutatorTerm<T, U> {
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

impl<T: Arith, U: Clone + Debug + PartialEq + PartialOrd + Ord + Hash> Mul<T> for CommutatorTerm<T, U> {
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

impl<U: Clone + Debug + PartialEq + PartialOrd + Ord + Hash> Mul<CommutatorTerm<NotNan<f32>, U>> for NotNan<f32> {
    type Output = CommutatorTerm<NotNan<f32>, U>;

    fn mul(self, rhs: CommutatorTerm<NotNan<f32>, U>) -> Self::Output {
        rhs.mul(self)
    }
}

impl<U: Clone + Debug + PartialEq + PartialOrd + Ord + Hash> Mul<CommutatorTerm<NotNan<f64>, U>> for NotNan<f64> {
    type Output = CommutatorTerm<NotNan<f64>, U>;

    fn mul(self, rhs: CommutatorTerm<NotNan<f64>, U>) -> Self::Output {
        rhs.mul(self)
    }
}

impl<U: Clone + Debug + PartialEq + PartialOrd + Ord + Hash> Mul<CommutatorTerm<i8, U>> for i8 {
    type Output = CommutatorTerm<i8, U>;

    fn mul(self, rhs: CommutatorTerm<i8, U>) -> Self::Output {
        rhs.mul(self)
    }
}

impl<U: Clone + Debug + PartialEq + PartialOrd + Ord + Hash> Mul<CommutatorTerm<i16, U>> for i16 {
    type Output = CommutatorTerm<i16, U>;

    fn mul(self, rhs: CommutatorTerm<i16, U>) -> Self::Output {
        rhs.mul(self)
    }
}

impl<U: Clone + Debug + PartialEq + PartialOrd + Ord + Hash> Mul<CommutatorTerm<i32, U>> for i32 {
    type Output = CommutatorTerm<i32, U>;

    fn mul(self, rhs: CommutatorTerm<i32, U>) -> Self::Output {
        rhs.mul(self)
    }
}

impl<U: Clone + Debug + PartialEq + PartialOrd + Ord + Hash> Mul<CommutatorTerm<i64, U>> for i64 {
    type Output = CommutatorTerm<i64, U>;

    fn mul(self, rhs: CommutatorTerm<i64, U>) -> Self::Output {
        rhs.mul(self)
    }
}

impl<U: Clone + Debug + PartialEq + PartialOrd + Ord + Hash> Mul<CommutatorTerm<isize, U>> for isize {
    type Output = CommutatorTerm<isize, U>;

    fn mul(self, rhs: CommutatorTerm<isize, U>) -> Self::Output {
        rhs.mul(self)
    }
}

impl<U: Clone + Debug + PartialEq + PartialOrd + Ord + Hash> Mul<CommutatorTerm<i128, U>> for i128 {
    type Output = CommutatorTerm<i128, U>;

    fn mul(self, rhs: CommutatorTerm<i128, U>) -> Self::Output {
        rhs.mul(self)
    }
}

impl<T: Arith, U: Clone + Debug + PartialEq + PartialOrd + Ord + Hash> Neg for CommutatorTerm<T, U> {
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


impl<T: Arith, U: Clone + Debug + PartialEq + PartialOrd + Ord + Hash> CommutatorTerm<T, U> {
    pub fn lyndon_sort(&mut self) {
        match self {
            // Do nothing, already sorted
            CommutatorTerm::Atom { .. } => {}
            CommutatorTerm::Expression{coefficient,  left, right} => {
                left.lyndon_sort();
                right.lyndon_sort();
                match left.cmp(&right) {
                    std::cmp::Ordering::Equal => *coefficient = T::zero(),
                    std::cmp::Ordering::Greater => {
                        *coefficient = -coefficient.clone();
                        std::mem::swap(left, right);
                    }
                    std::cmp::Ordering::Less => {}
                }
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

    /// Get the degree of the expression
    pub fn degree(&self) -> usize {
        match self {
            CommutatorTerm::Atom { .. } => 1,
            CommutatorTerm::Expression {  left, right , ..} => {
                left.degree() + right.degree()
            },
        }
    }

    pub fn is_zero(&self) -> bool {
        match self {
            CommutatorTerm::Atom { coefficient, ..} |
            CommutatorTerm::Expression { coefficient, .. } => coefficient.is_zero(),
        }
    }


    /// Returns the term with the unit coefficient
    #[must_use]
    pub fn unit(&self) -> Self {
        match self {
            Self::Atom {  atom, .. } => Self::Atom {
                coefficient: T::one(),
                atom: atom.clone()
            },
            Self::Expression { left, right, ..} =>
            Self::Expression {
                coefficient: T::one(),
                left: left.clone(),
                right: right.clone() }
        }
    }

    fn coefficient(&self) -> &T {
        match self {
            Self::Atom { coefficient, ..} | Self::Expression { coefficient, ..} => coefficient
        }
    }


    fn coefficient_mut(&mut self) -> &mut T {
        match self {
            Self::Atom { coefficient, ..} | Self::Expression { coefficient, ..} => coefficient
        }
    }

    fn left(&self) -> Option<&Self> {
        match self {
            CommutatorTerm::Atom { .. } => None,
            CommutatorTerm::Expression { left, .. } => Some(left),
        }
    }

    fn right(&self) -> Option<&Self> {
        match self {
            CommutatorTerm::Atom { .. } => None,
            CommutatorTerm::Expression { right, .. } => Some(right),
        }
    }

    fn left_mut(&mut self) -> Option<&mut Self> {
        match self {
            CommutatorTerm::Atom { .. } => None,
            CommutatorTerm::Expression { left, .. } => Some(left),
        }
    }

    fn right_mut(&mut self) -> Option<&mut Self> {
        match self {
            CommutatorTerm::Atom { .. } => None,
            CommutatorTerm::Expression { right, .. } => Some(right),
        }
    }

    fn find_decomposition_subterm_mut(&mut self, lyndon_basis_set: &HashSet<Self>) -> Option<&mut Self> {
        if let Self::Atom {..} = self {
            return None;
        }
        if lyndon_basis_set.contains(&self.unit()) {
            return None;
        }

        let Self::Expression { left, right, ..} = self else {
            return None;
        };

        if let Some(result) = left.find_decomposition_subterm_mut(lyndon_basis_set) {
            return Some(unsafe {
                std::ptr::from_mut(result).as_mut().unwrap()
            })
        }


        if let Some(result) = right.find_decomposition_subterm_mut(lyndon_basis_set) {
            return Some(unsafe {
                std::ptr::from_mut(result).as_mut().unwrap()
            })
        }

        Some(self)
    }

    #[must_use]
    pub fn lyndon_basis_decomposition(&self, lyndon_basis_set: &HashSet<Self>) -> Vec<Self> {
        if lyndon_basis_set.contains(&self.unit()) {
            return vec![self.clone()];
        }

        let mut lyndon_basis_terms = HashMap::<Self, Self>::new();
        let mut term_queue = vec![self.clone()];


        while let Some(mut t) = term_queue.pop() {
            t.lyndon_sort();
            if lyndon_basis_set.contains(&t.unit()) {
                lyndon_basis_terms.entry(t.unit()).and_modify(|x| *x.coefficient_mut() += t.coefficient().clone()).or_insert(t);
                continue;
            }
            let mut t1 = t.clone();
            let mut t2 = t.clone();
            let s1 = t1.find_decomposition_subterm_mut(lyndon_basis_set).unwrap();
            let s2 = t2.find_decomposition_subterm_mut(lyndon_basis_set).unwrap();


            if s1.is_zero() || s2.is_zero() || s1.left().unwrap() == s1.right().unwrap() {
                continue;
            }

            let (a, b) = {
                let s_prime = s1.left().unwrap();
                assert!(s_prime.left().is_some(), "Expected s' to be an expression, got: \n\ts' = {s_prime:?} \n\ts1 = {s1:?}");
                let a = s_prime.left().unwrap();
                let b = s_prime.right().unwrap();
                (a, b)
            };


            let s_dprime = s1.right().unwrap();

            let new_s_1 = Self::Expression {
                coefficient: s1.coefficient().clone(),
                left: Box::new(comm![a, s_dprime]),
                right: Box::new(b.clone())
            };

            let new_s_2 = Self::Expression {
                coefficient: s2.coefficient().clone(),
                left: Box::new(a.clone()),
                right: Box::new(comm![b, s_dprime])
            };

            *s1 = new_s_1;
            *s2 = new_s_2;

            if lyndon_basis_set.contains(&t1.unit()) {
                lyndon_basis_terms.entry(t1.unit()).and_modify(|x| *x.coefficient_mut() += t1.coefficient().clone()).or_insert(t1);
            } else {
                term_queue.push(t1);
            }

            if lyndon_basis_set.contains(&t2.unit()) {
                lyndon_basis_terms.entry(t2.unit()).and_modify(|x| *x.coefficient_mut() += t2.coefficient().clone()).or_insert(t2);
            } else {
                term_queue.push(t2);
            }
        }


        let mut lyndon_basis_terms = lyndon_basis_terms.into_values().collect::<Vec<_>>();

        lyndon_basis_terms.sort();
        lyndon_basis_terms
    }
}

impl<T: Arith, U: Clone + Debug + PartialEq + PartialOrd + Ord + Hash> Commutator<&Self>
    for CommutatorTerm<T, U>
{
    type Output = Self;

    fn commutator(&self, other: &Self) -> Self::Output {
        match (self, other) {
            (a @ Self::Atom {coefficient: c1, atom: a1}, 
                b @ Self::Atom {coefficient: c2, atom: a2}) => {
                let coefficient = if  a == b { T::zero() } else { c1.clone() * c2.clone() };
                let left = Box::new(Self::Atom { coefficient: T::one(), atom: a1.clone()});
                let right = Box::new(Self::Atom { coefficient: T::one(), atom: a2.clone()});
                Self::Expression{
                    coefficient,
                    left,
                    right,
                }
            }
            (Self::Atom {coefficient: c1, atom}, 
                Self::Expression {coefficient: c2, left: l1, right}) => {
                let coefficient = c1.clone() * c2.clone();
                let left = Box::new(Self::Atom { coefficient: T::one(), atom: atom.clone()});
                let right = Box::new(Self::Expression { coefficient: T::one(), left: l1.clone(), right: right.clone()});

                Self::Expression{
                    coefficient,
                    left,
                    right,
                }
            }
            (Self::Expression {coefficient: c1, left: l1, right: r1}, 
                Self::Atom{ coefficient: c2, atom }) => {
                let coefficient = c1.clone() * c2.clone();
                let left = Box::new(Self::Expression { coefficient: T::one(), left: l1.clone(), right: r1.clone() });
                let right = Box::new(Self::Atom { coefficient: T::one(), atom: atom.clone()});
                Self::Expression {
                    coefficient,
                    left,
                    right
                }
            }
            (a @ Self::Expression {coefficient: c1, left: l1, right: r1},
                b @ Self::Expression {coefficient: c2, left: l2, right: r2}) => {
                let coefficient = if a == b {
                    T::zero()
                } else {
                    c1.clone() * c2.clone()
                };
                let left = Box::new(Self::Expression { coefficient: T::one(), left: l1.clone(), right: r1.clone()});
                let right = Box::new(Self::Expression { coefficient: T::one(), left: l2.clone(), right: r2.clone()});

                Self::Expression {
                    coefficient,
                    left,
                    right,
                }
            }
        }
    }
}

impl<T: Arith, U: Clone + Debug + PartialEq + PartialOrd + Ord + Hash> PartialOrd
    for CommutatorTerm<T, U>
{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<T: Arith, U: Debug + Clone + Ord + Hash> Ord for CommutatorTerm<T, U> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match (self, other) {
            (Self::Atom { atom: a1, ..}, Self::Atom { atom: a2, ..}) => a1.cmp(a2),
            (Self::Atom { .. }, Self::Expression { left, .. }) => {
                match self.cmp(left) {
                    std::cmp::Ordering::Equal => std::cmp::Ordering::Less,
                    o => o,
                }
            }
            (Self::Expression { left, ..}, Self::Atom { .. }) => {
                match (**left).cmp(other) {
                    std::cmp::Ordering::Equal => std::cmp::Ordering::Greater,
                    o => o,
                }
            }
            (Self::Expression {left: l1, right: r1, ..}, 
                Self::Expression {left: l2, right: r2, ..}) => {
                match l1.cmp(l2) {
                    std::cmp::Ordering::Equal => r1.cmp(r2),
                    o => o,
                }
            }
        }
    }
}


impl<T: Arith, U: Generator + Debug + Clone + Ord + Hash> From<&LyndonWord<U>>
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

        let mut result = Self::Expression { coefficient: T::one(), left, right };
        result.lyndon_sort();
        *result.coefficient_mut() = T::one();
        result
    }
}


/// This can represent terms of the form `αe_1⋅⋅⋅e_n`
#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct FormalIndeterminate<T: Clone + Debug + Eq + Ord + PartialEq + PartialOrd + Hash, U: Arith> {
    coefficient: U,
    symbols: Vec<T>,
}

impl<T: Clone + Debug + Display + Eq + Ord + PartialEq + PartialOrd + Hash, U: Arith + Display> Display for FormalIndeterminate<T, U> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.coefficient)?;
        if !self.symbols.is_empty() {
            write!(f, " * ", )?;
            for symbol in &self.symbols {
                write!(f, "{symbol}")?;
            }
        }
        Ok(())
    }
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
    use crate::lyndon::{LyndonBasis, LyndonWordError};

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
        CommutatorTerm::from('A'),
        CommutatorTerm::from('B'),
        CommutatorTerm::Expression {
        coefficient: 1,
        left: Box::new(CommutatorTerm::from('A')),
        right: Box::new(CommutatorTerm::from('B')),
    })]
    #[case(
        CommutatorTerm::from('B'),
        CommutatorTerm::from('A'),
        CommutatorTerm::Expression {
            coefficient: 1,
            left: Box::new(CommutatorTerm::from('B')),
            right: Box::new(CommutatorTerm::from('A')),
        })]
    #[case(
        CommutatorTerm::Expression {
            coefficient: 2,
            left: Box::new(CommutatorTerm::from('A')),
            right: Box::new(CommutatorTerm::from('B'))
        },
        CommutatorTerm::Expression {
            coefficient: 3,
            left: Box::new(CommutatorTerm::from('B')),
            right: Box::new(CommutatorTerm::from('A'))
        },
        CommutatorTerm::Expression {
            coefficient: 6,
            left: Box::new(CommutatorTerm::Expression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::from('A')),
                right: Box::new(CommutatorTerm::from('B'))
            }),
            right: Box::new(CommutatorTerm::Expression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::from('B')),
                right: Box::new(CommutatorTerm::from('A'))
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
    #[case("AC", "ABB", Ordering::Less)]
    fn test_commutator_expression_ordering(
        #[case] word_1: &str,
        #[case] word_2: &str,
        #[case] expected_ordering: Ordering,
    ) -> Result<(), LyndonWordError> {
        let word_1 = word_1.parse::<LyndonWord<char>>()?;
        let word_2 = word_2.parse::<LyndonWord<char>>()?;
        let exp_1 = CommutatorTerm::<i32, char>::from(&word_1);
        let exp_2 = CommutatorTerm::from(&word_2);
        assert_eq!(exp_1.cmp(&exp_2), expected_ordering);

        Ok(())
    }

    #[rstest]
    #[case(CommutatorTerm::from('A'), CommutatorTerm::from('B'), Ordering::Less)]
    #[case(
        CommutatorTerm::from('B'),
        CommutatorTerm::from('A'),
        Ordering::Greater
    )]
    #[case(CommutatorTerm::from('A'), CommutatorTerm::from('A'), Ordering::Equal)]
    #[case(CommutatorTerm::from('A'), CommutatorTerm::from('A'), Ordering::Equal)]
    #[case(
        CommutatorTerm::from('A'),
        CommutatorTerm::Expression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::from('A')),
                right: Box::new(CommutatorTerm::from('B'))
            },
        Ordering::Less)]
    #[case(
        CommutatorTerm::from('B'),
        CommutatorTerm::Expression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::from('A')),
                right: Box::new(CommutatorTerm::from('B'))
            },
        Ordering::Greater)]
    #[case(
        CommutatorTerm::Expression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::from('A')),
                right: Box::new(CommutatorTerm::from('B'))
            },
        CommutatorTerm::from('A'),
        Ordering::Greater)]
    #[case(
        CommutatorTerm::Expression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::from('A')),
                right: Box::new(CommutatorTerm::from('B'))
            },
        CommutatorTerm::from('B'),
        Ordering::Less)]
    #[case(
        CommutatorTerm::Expression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::from('A')),
                right: Box::new(CommutatorTerm::from('B'))
            },
        CommutatorTerm::from('A'),
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
        comm![
            CommutatorTerm::from('B'),
            CommutatorTerm::from('A')
        ],
        -comm![
            CommutatorTerm::from('A'),
            CommutatorTerm::from('B')
        ]
    )]
    #[case(
        comm![
            CommutatorTerm::from('A'),
            CommutatorTerm::from('B')
        ],
        comm![
            CommutatorTerm::from('A'),
            CommutatorTerm::from('B')
        ],
    )]
    #[case(
        comm![
            comm![
                CommutatorTerm::from('B'),
                CommutatorTerm::from('A')
            ],
            CommutatorTerm::from('B')
        ],
        -comm![
            comm![
                CommutatorTerm::from('A'),
                CommutatorTerm::from('B')
            ],
            CommutatorTerm::from('B')
        ]
    )]
    #[case(
        comm![
            comm![
                CommutatorTerm::from('B'),
                CommutatorTerm::from('A')
            ],
            CommutatorTerm::from('A')
        ],
        comm![
            CommutatorTerm::from('A'),
            comm![
                CommutatorTerm::from('A'),
                CommutatorTerm::from('B')
            ]
        ],
    )]
    #[case(
        comm![CommutatorTerm::from('A'), CommutatorTerm::from('A')],
        0 * comm![CommutatorTerm::<i128, char>::from('A'), CommutatorTerm::from('A')]
    )]
    #[case(
        comm![
            comm![
                CommutatorTerm::from('B'),
                CommutatorTerm::from('A')
            ],
            comm![
                CommutatorTerm::from('A'),
                CommutatorTerm::from('B')
            ]
        ],
        0 * comm![
            comm![
                CommutatorTerm::<i128, char>::from('A'),
                CommutatorTerm::from('B')
            ],
            comm![
                CommutatorTerm::from('A'),
                CommutatorTerm::from('B')
            ]
        ],
    )]
    fn test_commutator_term_lyndon_sorting(
        #[case] mut term: CommutatorTerm<i128, char>,
        #[case] expected_term: CommutatorTerm<i128, char>,
    ) {
        term.lyndon_sort();
        println!("{term}");
        println!("{expected_term}");
        assert_eq!(term, expected_term);
    }

    #[rstest]
    #[case(
        comm![
            comm![
                CommutatorTerm::from('A'),
                CommutatorTerm::from('B')],
            CommutatorTerm::from('C')
        ],
        comm![
            CommutatorTerm::from('A'),
            comm![
                CommutatorTerm::from('B'),
                CommutatorTerm::from('C')
            ]
        ],
        comm![
            comm![
                CommutatorTerm::from('A'),
                CommutatorTerm::from('C')
            ],
        CommutatorTerm::from('B')]
    )]
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
        3,
        4,
        14 * comm![
            comm![
                comm![
                    CommutatorTerm::<i128, char>::from('A'),
                    CommutatorTerm::from('B')
                ],
                CommutatorTerm::from('B')
            ],
            CommutatorTerm::from('C')
        ],
        vec![
            14 * comm![
                comm![
                    comm![
                        CommutatorTerm::<i128, char>::from('A'),
                        CommutatorTerm::from('C')
                    ],
                    CommutatorTerm::from('B')
                ],
                CommutatorTerm::from('B')
            ],
            28 * comm![
                comm![
                    CommutatorTerm::<i128, char>::from('A'),
                    comm![
                        CommutatorTerm::from('B'),
                        CommutatorTerm::from('C')
                    ]
                ],
                CommutatorTerm::from('B')
            ],
            14 * comm![
                CommutatorTerm::<i128, char>::from('A'),
                comm![
                    CommutatorTerm::from('B'),
                    comm![
                        CommutatorTerm::from('B'),
                        CommutatorTerm::from('C')
                    ]
                ]
            ]
        ],
    )]
    #[case(
        4,
        5,
        - comm![
            comm![
                CommutatorTerm::from('A'),
                comm![
                    comm![
                        CommutatorTerm::from('A'),
                        CommutatorTerm::from('B')
                    ],
                    CommutatorTerm::from('B')
                ]
            ],
            CommutatorTerm::from('C')
        ],
        vec![]
    )]
    fn test_commutator_decomposition(
        #[case] num_generators: usize,
        #[case] max_degree: usize,
        #[case] term: CommutatorTerm<i128, char>,
        #[case] mut expected_basis_terms: Vec<CommutatorTerm<i128, char>>
    ) {
            use crate::lyndon::Sort;

        expected_basis_terms.sort();
        println!("term: {term}");
        let basis = LyndonBasis::<char>::new(num_generators, Sort::Lexicographical);
        let basis_set = basis
            .generate_basis(max_degree)
            .iter()
            .map(CommutatorTerm::<i128, char>::from)
            .collect::<HashSet<_>>();
        

        let mut basis_vec = basis_set.iter().collect::<Vec<_>>();
        basis_vec.sort();
        println!("Commutator Basis");
        for (i, term) in basis_vec.into_iter().enumerate() {
            let basis_str = format!("{term}");
            if term.degree() == 5 && !basis_str.contains('D') {
                println!("{i}: {term}");
            }
        }
        println!();

        let basis_terms = term.lyndon_basis_decomposition(&basis_set);
        assert_eq!(basis_terms.len(), expected_basis_terms.len());
        for (basis_term, expected_basis_term) in basis_terms.iter().zip(&expected_basis_terms) {
            assert_eq!(basis_term, expected_basis_term, "{basis_term} != {expected_basis_term}");
        }
    }

    #[rstest]
    #[case(
        CommutatorTerm::from('A'),
        CommutatorTerm::from('B'),
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
        comm![CommutatorTerm::from('A'), CommutatorTerm::from('B')],
        comm![CommutatorTerm::from('C'), CommutatorTerm::Atom { atom: 'D', coefficient: 1 }],
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
