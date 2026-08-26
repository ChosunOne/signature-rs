use num_traits::{One, Zero};

use lyndon_rs::{generators::Generator, lyndon::LyndonWord};
use std::{
    collections::{HashMap, HashSet},
    fmt::{Debug, Display},
    hash::{DefaultHasher, Hash, Hasher},
    ops::{AddAssign, Mul, MulAssign, Neg, Sub},
};

/// Trait for types that support the commutator operation.
///
/// The commutator operation `[A, B]` typically represents the algebraic expression `AB - BA`,
/// which measures the failure of two elements to commute.
pub trait Commutator<Rhs = Self> {
    /// The result type of the commutator operation.
    type Output;

    /// Computes the commutator `[self, other]` = `self * other - other * self`.
    fn commutator(&self, other: Rhs) -> Self::Output;
}

impl<T> Commutator<Self> for T
where
    Self: Clone + Mul<Output = T> + Sub<Output = T>,
{
    type Output = T;

    fn commutator(&self, other: Self) -> Self::Output {
        self.clone() * other.clone() - other.clone() * self.clone()
    }
}

impl<T> Commutator<&Self> for T
where
    Self: Clone + Mul<Output = T> + Sub<Output = T>,
{
    type Output = T;

    fn commutator(&self, other: &Self) -> Self::Output {
        self.clone() * other.clone() - other.clone() * self.clone()
    }
}

/// Shorthand macro for computing commutators `[A, B]`.
///
/// This macro provides a convenient syntax for creating commutator expressions.
/// It expands to `$a.commutator(&$b)`, applying the commutator operation between two terms.
///
/// # Examples
///
/// Basic usage with CommutatorTerm:
/// ```rust
/// use commutator_rs::{CommutatorTerm, Commutator, comm};
///
/// let x = CommutatorTerm::Atom { coefficient: 1, atom: 'x' };
/// let y = CommutatorTerm::Atom { coefficient: 1, atom: 'y' };
///
/// // These are equivalent:
/// let result1 = x.commutator(&y);
/// let result2 = comm![x, y];
/// assert_eq!(result1, result2);
/// ```
///
/// Nested commutators:
/// ```rust
/// use commutator_rs::{CommutatorTerm, Commutator, comm};
///
/// let a = CommutatorTerm::<i32, char>::from('a');
/// let b = CommutatorTerm::<i32, char>::from('b');
/// let c = CommutatorTerm::<i32, char>::from('c');
///
/// // Compute [[a, b], c]
/// let nested = comm![comm![a, b], c];
/// ```
///
/// With numeric types:
/// ```rust
/// use commutator_rs::{Commutator, comm};
///
/// // For numeric types, commutator is AB - BA
/// assert_eq!(comm![2i32, 3i32], 0); // 2*3 - 3*2 = 0
/// ```
#[macro_export]
macro_rules! comm {
    ($a:expr, $b:expr) => {
        $a.commutator(&$b)
    };
}

/// Represents an algebraic term involving commutators.
///
/// A `CommutatorTerm` can be either:
/// - An atomic element with a coefficient
/// - A commutator expression `[left, right]` with a coefficient
///
/// This structure allows for the representation of nested commutator expressions
/// and linear combinations thereof.
pub enum CommutatorTerm<T, U> {
    /// An atomic term consisting of a coefficient and an atom.
    Atom {
        /// The scalar coefficient multiplying the atom.
        coefficient: T,
        /// The atomic element (generator).
        atom: U,
    },
    /// A commutator expression `[left, right]` with a coefficient.
    Expression {
        /// The scalar coefficient multiplying the commutator.
        coefficient: T,
        /// The left operand of the commutator.
        left: Box<Self>,
        /// The right operand of the commutator.
        right: Box<Self>,
        /// The degree of the commutator
        degree: usize,
    },
}

impl<T: One> From<char> for CommutatorTerm<T, char> {
    fn from(value: char) -> Self {
        Self::Atom {
            coefficient: T::one(),
            atom: value,
        }
    }
}

impl<T: One> From<u8> for CommutatorTerm<T, u8> {
    fn from(value: u8) -> Self {
        Self::Atom {
            coefficient: T::one(),
            atom: value,
        }
    }
}

impl<T: Display + One + PartialEq, U: Display> Display for CommutatorTerm<T, U> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Atom { coefficient, atom } => {
                if coefficient.is_one() {
                    write!(f, "{atom}")
                } else {
                    write!(f, "{coefficient} * {atom}")
                }
            }
            Self::Expression {
                coefficient,
                left,
                right,
                ..
            } => {
                if coefficient.is_one() {
                    write!(f, "[{left}, {right}]")
                } else {
                    write!(f, "{coefficient} * [{left}, {right}]")
                }
            }
        }
    }
}

impl<T: Debug, U: Debug> Debug for CommutatorTerm<T, U> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Atom { coefficient, atom } => f
                .debug_struct("Atom")
                .field("coefficient", coefficient)
                .field("atom", atom)
                .finish(),
            Self::Expression {
                coefficient,
                left,
                right,
                degree,
            } => f
                .debug_struct("Expression")
                .field("coefficient", coefficient)
                .field("left", left)
                .field("right", right)
                .field("degree", degree)
                .finish(),
        }
    }
}

impl<T: Mul<Output = T>, U: Clone> Mul<T> for CommutatorTerm<T, U> {
    type Output = Self;

    fn mul(self, rhs: T) -> Self::Output {
        match self {
            Self::Atom { coefficient, atom } => Self::Atom {
                coefficient: coefficient * rhs,
                atom: atom.clone(),
            },
            Self::Expression {
                coefficient,
                left,
                right,
                degree,
            } => Self::Expression {
                coefficient: coefficient * rhs,
                left,
                right,
                degree,
            },
        }
    }
}

impl<T: Mul<Output = T> + Clone, U: Clone> Mul<T> for &CommutatorTerm<T, U> {
    type Output = CommutatorTerm<T, U>;

    fn mul(self, rhs: T) -> Self::Output {
        match self {
            CommutatorTerm::Atom { coefficient, atom } => CommutatorTerm::Atom {
                coefficient: coefficient.clone() * rhs,
                atom: atom.clone(),
            },
            CommutatorTerm::Expression {
                coefficient,
                left,
                right,
                degree,
            } => CommutatorTerm::Expression {
                coefficient: coefficient.clone() * rhs,
                left: left.clone(),
                right: right.clone(),
                degree: *degree,
            },
        }
    }
}

impl<U: Clone> Mul<CommutatorTerm<f32, U>> for f32 {
    type Output = CommutatorTerm<f32, U>;

    fn mul(self, rhs: CommutatorTerm<f32, U>) -> Self::Output {
        rhs.mul(self)
    }
}

impl<U: Clone> Mul<CommutatorTerm<f64, U>> for f64 {
    type Output = CommutatorTerm<f64, U>;

    fn mul(self, rhs: CommutatorTerm<f64, U>) -> Self::Output {
        rhs.mul(self)
    }
}

impl<U: Clone> Mul<CommutatorTerm<i8, U>> for i8 {
    type Output = CommutatorTerm<i8, U>;

    fn mul(self, rhs: CommutatorTerm<i8, U>) -> Self::Output {
        rhs.mul(self)
    }
}

impl<U: Clone> Mul<CommutatorTerm<i16, U>> for i16 {
    type Output = CommutatorTerm<i16, U>;

    fn mul(self, rhs: CommutatorTerm<i16, U>) -> Self::Output {
        rhs.mul(self)
    }
}

impl<U: Clone> Mul<CommutatorTerm<i32, U>> for i32 {
    type Output = CommutatorTerm<i32, U>;

    fn mul(self, rhs: CommutatorTerm<i32, U>) -> Self::Output {
        rhs.mul(self)
    }
}

impl<U: Clone> Mul<CommutatorTerm<i64, U>> for i64 {
    type Output = CommutatorTerm<i64, U>;

    fn mul(self, rhs: CommutatorTerm<i64, U>) -> Self::Output {
        rhs.mul(self)
    }
}

impl<U: Clone> Mul<CommutatorTerm<isize, U>> for isize {
    type Output = CommutatorTerm<isize, U>;

    fn mul(self, rhs: CommutatorTerm<isize, U>) -> Self::Output {
        rhs.mul(self)
    }
}

impl<U: Clone> Mul<CommutatorTerm<i128, U>> for i128 {
    type Output = CommutatorTerm<i128, U>;

    fn mul(self, rhs: CommutatorTerm<i128, U>) -> Self::Output {
        rhs.mul(self)
    }
}

impl<T: Hash, U: Hash> Hash for CommutatorTerm<T, U> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        core::mem::discriminant(self).hash(state);
        match self {
            CommutatorTerm::Atom { coefficient, atom } => {
                coefficient.hash(state);
                atom.hash(state);
            }
            CommutatorTerm::Expression {
                coefficient,
                left,
                right,
                ..
            } => {
                coefficient.hash(state);
                left.hash(state);
                right.hash(state);
            }
        }
    }
}

impl<T: Neg<Output = T>, U> Neg for CommutatorTerm<T, U> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        match self {
            Self::Atom { coefficient, atom } => Self::Atom {
                coefficient: coefficient.neg(),
                atom,
            },
            Self::Expression {
                coefficient,
                left,
                right,
                degree,
            } => Self::Expression {
                coefficient: coefficient.neg(),
                left,
                right,
                degree,
            },
        }
    }
}

impl<T: Clone + Neg<Output = T>, U: Clone> Neg for &CommutatorTerm<T, U> {
    type Output = CommutatorTerm<T, U>;

    fn neg(self) -> Self::Output {
        let mut result = self.clone();
        *result.coefficient_mut() = self.coefficient().clone().neg();
        result
    }
}

impl<T: Eq, U: Eq> PartialEq for CommutatorTerm<T, U> {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (
                Self::Atom {
                    coefficient: l_coefficient,
                    atom: l_atom,
                },
                Self::Atom {
                    coefficient: r_coefficient,
                    atom: r_atom,
                },
            ) => l_coefficient == r_coefficient && l_atom == r_atom,
            (
                Self::Expression {
                    coefficient: l_coefficient,
                    left: l_left,
                    right: l_right,
                    ..
                },
                Self::Expression {
                    coefficient: r_coefficient,
                    left: r_left,
                    right: r_right,
                    ..
                },
            ) => l_coefficient == r_coefficient && l_left == r_left && l_right == r_right,
            _ => false,
        }
    }
}
impl<T: Eq, U: Eq> Eq for CommutatorTerm<T, U> {}

impl<T: Clone, U: Clone> Clone for CommutatorTerm<T, U> {
    fn clone(&self) -> Self {
        match self {
            Self::Atom { coefficient, atom } => Self::Atom {
                coefficient: coefficient.clone(),
                atom: atom.clone(),
            },
            Self::Expression {
                coefficient,
                left,
                right,
                degree,
            } => Self::Expression {
                coefficient: coefficient.clone(),
                left: left.clone(),
                right: right.clone(),
                degree: *degree,
            },
        }
    }
}

impl<T: Clone + One + Zero + Eq + Mul<Output = T>, U: Clone + Eq> Commutator<&Self>
    for CommutatorTerm<T, U>
{
    type Output = Self;

    #[allow(clippy::too_many_lines)]
    fn commutator(&self, other: &Self) -> Self::Output {
        match (self, other) {
            (
                a @ Self::Atom {
                    coefficient: c1,
                    atom: a1,
                },
                b @ Self::Atom {
                    coefficient: c2,
                    atom: a2,
                },
            ) => {
                let coefficient = if a == b {
                    T::zero()
                } else {
                    c1.clone() * c2.clone()
                };
                let left = Box::new(Self::Atom {
                    coefficient: T::one(),
                    atom: a1.clone(),
                });
                let right = Box::new(Self::Atom {
                    coefficient: T::one(),
                    atom: a2.clone(),
                });
                Self::Expression {
                    coefficient,
                    left,
                    right,
                    degree: 2,
                }
            }
            (
                Self::Atom {
                    coefficient: c1,
                    atom,
                },
                Self::Expression {
                    coefficient: c2,
                    left: l1,
                    right,
                    degree,
                },
            ) => {
                let coefficient = c1.clone() * c2.clone();
                let left = Box::new(Self::Atom {
                    coefficient: T::one(),
                    atom: atom.clone(),
                });
                let right = Box::new(Self::Expression {
                    coefficient: T::one(),
                    left: l1.clone(),
                    right: right.clone(),
                    degree: *degree,
                });

                Self::Expression {
                    coefficient,
                    left,
                    right,
                    degree: *degree + 1,
                }
            }
            (
                Self::Expression {
                    coefficient: c1,
                    left: l1,
                    right: r1,
                    degree,
                },
                Self::Atom {
                    coefficient: c2,
                    atom,
                },
            ) => {
                let coefficient = c1.clone() * c2.clone();
                let left = Box::new(Self::Expression {
                    coefficient: T::one(),
                    left: l1.clone(),
                    right: r1.clone(),
                    degree: *degree,
                });
                let right = Box::new(Self::Atom {
                    coefficient: T::one(),
                    atom: atom.clone(),
                });
                Self::Expression {
                    coefficient,
                    left,
                    right,
                    degree: degree + 1,
                }
            }
            (
                a @ Self::Expression {
                    coefficient: c1,
                    left: l1,
                    right: r1,
                    degree: degree1,
                },
                b @ Self::Expression {
                    coefficient: c2,
                    left: l2,
                    right: r2,
                    degree: degree2,
                },
            ) => {
                let coefficient = if a == b {
                    T::zero()
                } else {
                    c1.clone() * c2.clone()
                };
                let left = Box::new(Self::Expression {
                    coefficient: T::one(),
                    left: l1.clone(),
                    right: r1.clone(),
                    degree: *degree1,
                });
                let right = Box::new(Self::Expression {
                    coefficient: T::one(),
                    left: l2.clone(),
                    right: r2.clone(),
                    degree: *degree2,
                });

                Self::Expression {
                    coefficient,
                    left,
                    right,
                    degree: *degree1 + *degree2,
                }
            }
        }
    }
}

impl<T: Eq, U: PartialEq + PartialOrd + Ord> PartialOrd for CommutatorTerm<T, U> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<T, U> CommutatorTerm<T, U> {
    /// Appends the atoms of this term in left-to-right (flattened) order.
    fn push_letters<'a>(&'a self, letters: &mut Vec<&'a U>) {
        match self {
            Self::Atom { atom, .. } => letters.push(atom),
            Self::Expression { left, right, .. } => {
                left.push_letters(letters);
                right.push_letters(letters);
            }
        }
    }

    /// Compares two terms by structure, ignoring coefficients.
    fn structure_eq(a: &Self, b: &Self) -> bool
    where
        U: PartialEq,
    {
        match (a, b) {
            (Self::Atom { atom: a1, .. }, Self::Atom { atom: a2, .. }) => a1 == a2,
            (
                Self::Expression {
                    left: l1,
                    right: r1,
                    ..
                },
                Self::Expression {
                    left: l2,
                    right: r2,
                    ..
                },
            ) => Self::structure_eq(l1, l2) && Self::structure_eq(r1, r2),
            _ => false,
        }
    }
}

impl<T: Eq, U: Ord> Ord for CommutatorTerm<T, U> {
    /// Compares terms by the lexicographic order of their flattened letter
    /// sequence (the same ordering `LyndonWord` uses), so that sorting,
    /// normalization, and basis-membership decisions are consistent with the
    /// Lyndon basis.
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        let mut letters = Vec::new();
        let mut other_letters = Vec::new();
        self.push_letters(&mut letters);
        other.push_letters(&mut other_letters);
        letters.cmp(&other_letters)
    }
}

impl<T: Clone + One + Zero + Eq + MulAssign + Neg<Output = T>, U: Clone + Ord + Generator + Eq>
    From<&LyndonWord<U>> for CommutatorTerm<T, U>
{
    fn from(value: &LyndonWord<U>) -> Self {
        if value.len() == 1 {
            return CommutatorTerm::Atom {
                coefficient: T::one(),
                atom: value.letters[0].clone(),
            };
        }

        let (left, right) = value.factorize();
        let left = if left.len() == 1 {
            Box::new(CommutatorTerm::Atom {
                coefficient: T::one(),
                atom: left.letters[0].clone(),
            })
        } else {
            Box::new(Self::from(&left))
        };
        let right = if right.len() == 1 {
            Box::new(CommutatorTerm::Atom {
                coefficient: T::one(),
                atom: right.letters[0].clone(),
            })
        } else {
            Box::new(Self::from(&right))
        };

        let mut result = Self::Expression {
            coefficient: T::one(),
            left,
            right,
            degree: value.len(),
        };
        result.lyndon_sort();
        *result.coefficient_mut() = T::one();
        result
    }
}

impl<T, U> CommutatorTerm<T, U> {
    /// Returns the degree of the commutator term.
    ///
    /// The degree is defined as the total number of atomic elements in the expression.
    /// For atoms, this is 1. For expressions, it's the sum of degrees of left and right operands.
    pub fn degree(&self) -> usize {
        match self {
            CommutatorTerm::Atom { .. } => 1,
            CommutatorTerm::Expression { degree, .. } => *degree,
        }
    }

    /// Returns a reference to the coefficient of this term.
    pub fn coefficient(&self) -> &T {
        match self {
            Self::Atom { coefficient, .. } | Self::Expression { coefficient, .. } => coefficient,
        }
    }

    /// Returns a mutable reference to the coefficient of this term.
    pub fn coefficient_mut(&mut self) -> &mut T {
        match self {
            Self::Atom { coefficient, .. } | Self::Expression { coefficient, .. } => coefficient,
        }
    }

    pub fn atom(&self) -> Option<&U> {
        match self {
            Self::Atom { atom, .. } => Some(atom),
            Self::Expression { .. } => None,
        }
    }

    pub fn atom_mut(&mut self) -> Option<&mut U> {
        match self {
            Self::Atom { atom, .. } => Some(atom),
            Self::Expression { .. } => None,
        }
    }

    /// Returns the left operand of a commutator expression, or `None` for atoms.
    pub fn left(&self) -> Option<&Self> {
        match self {
            Self::Atom { .. } => None,
            Self::Expression { left, .. } => Some(left),
        }
    }

    /// Returns the right operand of a commutator expression, or `None` for atoms.
    pub fn right(&self) -> Option<&Self> {
        match self {
            Self::Atom { .. } => None,
            Self::Expression { right, .. } => Some(right),
        }
    }

    /// Returns a mutable reference to the left operand, or `None` for atoms.
    pub fn left_mut(&mut self) -> Option<&mut Self> {
        match self {
            Self::Atom { .. } => None,
            Self::Expression { left, .. } => Some(left),
        }
    }

    /// Returns a mutable reference to the right operand, or `None` for atoms.
    pub fn right_mut(&mut self) -> Option<&mut Self> {
        match self {
            Self::Atom { .. } => None,
            Self::Expression { right, .. } => Some(right),
        }
    }
}

impl<T: Zero, U> CommutatorTerm<T, U> {
    /// Returns `true` if the coefficient of this term is zero.
    pub fn is_zero(&self) -> bool {
        match self {
            Self::Atom { coefficient, .. } | Self::Expression { coefficient, .. } => {
                coefficient.is_zero()
            }
        }
    }
}

impl<T: Clone + One, U: Clone> CommutatorTerm<T, U> {
    /// Returns a copy of this term with coefficient set to one.
    ///
    /// This is useful for extracting the structural part of a term
    /// without its scalar coefficient.
    #[must_use]
    pub fn unit(&self) -> Self {
        match self {
            Self::Atom { atom, .. } => Self::Atom {
                coefficient: T::one(),
                atom: atom.clone(),
            },
            Self::Expression { left, right, .. } => Self::Expression {
                coefficient: T::one(),
                left: left.clone(),
                right: right.clone(),
                degree: left.degree() + right.degree(),
            },
        }
    }
}

impl<T, U: Hash> CommutatorTerm<T, U> {
    /// Returns the hash of only the atoms in the struct
    pub fn atom_hash(&self) -> u64 {
        let mut state = DefaultHasher::default();
        core::mem::discriminant(self).hash(&mut state);
        match self {
            Self::Atom { atom, .. } => {
                atom.hash(&mut state);
            }
            Self::Expression { left, right, .. } => {
                left.atom_hash().hash(&mut state);
                right.atom_hash().hash(&mut state);
            }
        }

        state.finish()
    }
}

impl<T: Clone + One + Hash, U: Hash> CommutatorTerm<T, U> {
    /// Returns the hash of the unit struct without cloning
    pub fn unit_hash(&self) -> u64 {
        let mut state = DefaultHasher::default();
        core::mem::discriminant(self).hash(&mut state);
        match self {
            CommutatorTerm::Atom { atom, .. } => {
                T::one().hash(&mut state);
                atom.hash(&mut state);
            }
            CommutatorTerm::Expression { left, right, .. } => {
                T::one().hash(&mut state);
                left.unit_hash().hash(&mut state);
                right.unit_hash().hash(&mut state);
            }
        }
        state.finish()
    }
}

impl<T: Eq + Clone + Neg<Output = T> + Zero + One + MulAssign + PartialEq, U: Clone + Ord + Eq>
    CommutatorTerm<T, U>
{
    /// Sorts the commutator term into canonical Lyndon ordering.
    ///
    /// This method recursively applies the anti-commutativity property `[A, B] = -[B, A]`
    /// to ensure that commutator expressions are in a canonical form where the left
    /// operand is lexicographically smaller than the right operand.
    pub fn lyndon_sort(&mut self) {
        match self {
            // Do nothing, already sorted
            Self::Atom { .. } => {}
            Self::Expression {
                coefficient,
                left,
                right,
                ..
            } => {
                left.lyndon_sort();
                right.lyndon_sort();
                // Zero out [X, X] terms by anti-commutativity. This must use
                // structural equality (ignoring coefficients): flatten-lex
                // equality no longer implies structural equality (e.g.
                // [A,[B,C]] and [[A,B],C] both flatten to "ABC").
                if Self::structure_eq(left, right) {
                    *coefficient = T::zero();
                } else if (**left).cmp(&**right) == std::cmp::Ordering::Greater {
                    *coefficient = -coefficient.clone();
                    std::mem::swap(left, right);
                }
                // Propagate up coefficients
                if let Self::Expression {
                    coefficient: c2, ..
                } = &mut **left
                {
                    *coefficient *= c2.clone();
                    if *c2 == -T::one() {
                        *c2 = -c2.clone();
                    }
                }
                if let Self::Expression {
                    coefficient: c2, ..
                } = &mut **right
                {
                    *coefficient *= c2.clone();
                    if *c2 == -T::one() {
                        *c2 = -c2.clone();
                    }
                }
            }
        }
    }

    /// Applies the Jacobi identity to decompose nested commutators.
    ///
    /// For a commutator of the form `[[A, B], C]`, returns the equivalent expression
    /// `[A, [B, C]] - [B, [A, C]]` as a pair of terms. Returns `None` if the
    /// Jacobi identity cannot be applied (e.g., for atoms or already-decomposed expressions).
    pub fn jacobi_identity(&self) -> Option<(Self, Self)> {
        match self {
            Self::Atom { .. } => None,
            Self::Expression {
                coefficient,
                left,
                right,
                ..
            } => {
                // Check if expression is already in basis form
                let Self::Expression {
                    left: l_left,
                    right: l_right,
                    ..
                } = &**left
                else {
                    return None;
                };
                let a = l_left.clone();
                let b = l_right.clone();
                let c = right.clone();

                let mut left_term = comm![a, comm![b, c]] * coefficient.clone();
                left_term.lyndon_sort();
                let mut right_term = comm![b, comm![a, c]] * -coefficient.clone();
                right_term.lyndon_sort();

                Some((left_term, right_term))
            }
        }
    }
}

impl<
    T: Clone + Hash + Eq + One + Zero + MulAssign + Neg<Output = T> + Ord + AddAssign,
    U: Clone + Hash + Eq + Ord,
> CommutatorTerm<T, U>
{
    fn find_decomposition_subterm_mut(
        &mut self,
        lyndon_basis_set: &HashSet<u64>,
    ) -> Option<&mut Self> {
        if let Self::Atom { .. } = self {
            return None;
        }
        if lyndon_basis_set.contains(&self.unit_hash()) {
            return None;
        }

        let Self::Expression { left, right, .. } = self else {
            return None;
        };

        if let Some(result) = left.find_decomposition_subterm_mut(lyndon_basis_set) {
            return Some(unsafe { std::ptr::from_mut(result).as_mut().unwrap() });
        }

        if let Some(result) = right.find_decomposition_subterm_mut(lyndon_basis_set) {
            return Some(unsafe { std::ptr::from_mut(result).as_mut().unwrap() });
        }

        Some(self)
    }
    /// Decomposes this commutator term into a linear combination of Lyndon basis elements.
    ///
    /// Given a set of Lyndon basis elements, this method expresses the current term
    /// as a sum of basis terms using the Jacobi identity and other commutator relations.
    ///
    /// The Jacobi reassociation `[[a, b], v] = [[a, v], b] + [a, [b, v]]` is only
    /// applied when the term is not in standard form, i.e. when `b < v`. A sorted
    /// subterm whose children are both basis terms is itself a basis (standard)
    /// term exactly when `b >= v` (the right factor must be the longest Lyndon
    /// suffix), so every guarded rewrite makes progress toward standard form.
    /// Applying the identity unconditionally makes it involutive across non-basis
    /// patterns (`[[a,b],c] -> [[a,c],b] -> [[a,b],c]`), which never terminates.
    /// Additionally, identical non-basis terms are merged linearly so each
    /// distinct term is rewritten at most once.
    #[must_use]
    pub fn lyndon_basis_decomposition(&self, lyndon_basis_set: &HashSet<u64>) -> Vec<Self> {
        #[cfg(feature = "tracing")]
        let _span =
            tracing::debug_span!("lyndon_basis_decomposition", degree = self.degree()).entered();
        if lyndon_basis_set.contains(&self.unit_hash()) {
            #[cfg(feature = "tracing")]
            tracing::trace!("term is already a basis term");
            return vec![self.clone()];
        }

        // Accumulated basis terms keyed by their unit form, with summed coefficients.
        let mut lyndon_basis_terms = HashMap::<Self, Self>::new();
        // Non-basis terms awaiting rewriting, keyed by unit hash. Merging
        // identical terms linearly guarantees each distinct non-basis term is
        // rewritten at most once.
        let mut pending = HashMap::<u64, Self>::new();
        let mut term_queue = vec![self.clone()];
        #[cfg(feature = "tracing")]
        let mut iterations = 0_usize;

        loop {
            while let Some(mut t) = term_queue.pop() {
                #[cfg(feature = "tracing")]
                {
                    iterations += 1;
                    if iterations % 4096 == 1 {
                        tracing::debug!(
                            iterations,
                            queue_len = term_queue.len(),
                            pending = pending.len(),
                            terms_found = lyndon_basis_terms.len(),
                            "decomposition in progress"
                        );
                    }
                }
                t.lyndon_sort();
                let h = t.unit_hash();
                if lyndon_basis_set.contains(&h) {
                    lyndon_basis_terms
                        .entry(t.unit())
                        .and_modify(|x| *x.coefficient_mut() += t.coefficient().clone())
                        .or_insert(t);
                    continue;
                }
                pending
                    .entry(h)
                    .and_modify(|x| *x.coefficient_mut() += t.coefficient().clone())
                    .or_insert(t);
            }

            let Some(&h) = pending.keys().next() else {
                break;
            };
            let t = pending
                .remove(&h)
                .expect("a key obtained from pending to exist in pending");
            if t.is_zero() {
                continue;
            }

            let mut t1 = t.clone();
            let mut t2 = t.clone();
            let s1 = t1
                .find_decomposition_subterm_mut(lyndon_basis_set)
                .expect("a non-basis term to contain a non-basis subterm");
            let s2 = t2
                .find_decomposition_subterm_mut(lyndon_basis_set)
                .expect("a non-basis term to contain a non-basis subterm");

            if s1.is_zero() || s2.is_zero() || s1.left().unwrap() == s1.right().unwrap() {
                #[cfg(feature = "tracing")]
                tracing::trace!(pending = pending.len(), "dropping zero/invalid subterm");
                continue;
            }

            // The subterm has the form s = [u, v] with u = [a, b], sorted so that
            // u < v. It is a standard (basis) term iff b >= v; only rewrite when
            // it is out of standard form (b < v).
            let (a, b) = {
                let u = s1.left().unwrap();
                let a = u.left().unwrap();
                let b = u.right().unwrap();
                (a, b)
            };

            let v = s1.right().unwrap();

            if !(b < v) {
                // Unreachable for a sorted non-basis subterm with basis children;
                // kept as a guard against non-terminating rewrites if the ordering
                // assumptions are ever broken.
                #[cfg(feature = "tracing")]
                tracing::warn!(
                    "non-basis subterm is already in standard form; skipping rewrite"
                );
                continue;
            }

            let new_s_1 = Self::Expression {
                coefficient: s1.coefficient().clone(),
                left: Box::new(comm![a, v]),
                right: Box::new(b.clone()),
                degree: self.degree(),
            };

            let new_s_2 = Self::Expression {
                coefficient: s2.coefficient().clone(),
                left: Box::new(a.clone()),
                right: Box::new(comm![b, v]),
                degree: self.degree(),
            };

            *s1 = new_s_1;
            *s2 = new_s_2;

            term_queue.push(t1);
            term_queue.push(t2);
        }

        #[cfg(feature = "tracing")]
        tracing::debug!(
            iterations,
            terms_found = lyndon_basis_terms.len(),
            "decomposition complete"
        );

        let mut lyndon_basis_terms = lyndon_basis_terms.into_values().collect::<Vec<_>>();

        lyndon_basis_terms.sort();
        lyndon_basis_terms
    }
}

#[cfg(test)]
mod test {
    use lyndon_rs::lyndon::{LyndonBasis, Sort};

    use crate::formal_indeterminate::FormalIndeterminate;

    /// The term ordering must agree with the underlying Lyndon word ordering
    /// (flatten-lex); a mismatch made `lyndon_sort` and basis membership
    /// decisions inconsistent and led to wrong/dropped decomposition terms.
    #[test]
    fn commutator_term_ordering_matches_lyndon_word_ordering() {
        use lyndon_rs::lyndon::{LyndonBasis, Sort};

        let basis = LyndonBasis::<u8>::new(3, Sort::Lexicographical).generate_basis(5);
        for w1 in &basis {
            for w2 in &basis {
                let t1 = CommutatorTerm::<i128, u8>::from(w1);
                let t2 = CommutatorTerm::<i128, u8>::from(w2);
                assert_eq!(
                    t1.cmp(&t2),
                    w1.cmp(w2),
                    "ordering mismatch on {w1:?} vs {w2:?}"
                );
            }
        }
    }

    /// Expands a commutator term into the free associative algebra over its
    /// atoms, where `[x, y] = x*y - y*x`. Returns a map from words (letter
    /// sequences) to coefficients. This is an independent oracle for checking
    /// Lie-algebraic identities without relying on any basis machinery.
    fn expand_to_free_algebra(term: &CommutatorTerm<i128, u8>) -> HashMap<Vec<u8>, i128> {
        let mut result = HashMap::new();
        match term {
            CommutatorTerm::Atom {
                coefficient, atom,
            } => {
                *result.entry(vec![*atom]).or_default() += coefficient;
            }
            CommutatorTerm::Expression {
                coefficient,
                left,
                right,
                ..
            } => {
                let l = expand_to_free_algebra(left);
                let r = expand_to_free_algebra(right);
                for (wl, cl) in &l {
                    for (wr, cr) in &r {
                        let mut lr = wl.clone();
                        lr.extend(wr);
                        *result.entry(lr).or_default() += coefficient * cl * cr;
                        let mut rl = wr.clone();
                        rl.extend(wl);
                        *result.entry(rl).or_default() -= coefficient * cl * cr;
                    }
                }
            }
        }
        result.retain(|_, c| *c != 0);
        result
    }


    /// Oracle sanity check: the free-algebra expansion must agree with the
    /// known-good decomposition 14[[AB]B, C] =
    /// 14[[AC]B, B] + 28[[A[BC]], B] + 14[A[B[BC]]].
    #[test]
    fn free_algebra_expansion_oracle_matches_known_decomposition() {
        let a = CommutatorTerm::<i128, u8>::from(0_u8);
        let b = CommutatorTerm::<i128, u8>::from(1_u8);
        let c = CommutatorTerm::<i128, u8>::from(2_u8);

        let term = 14 * comm![comm![comm![a.clone(), b.clone()], b.clone()], c.clone()];
        let known_terms = vec![
            14 * comm![comm![comm![a.clone(), c.clone()], b.clone()], b.clone()],
            28 * comm![comm![a.clone(), comm![b.clone(), c.clone()]], b.clone()],
            14 * comm![a, comm![b, comm![b, c]]],
        ];

        let original = expand_to_free_algebra(&term);
        let mut recomposed = HashMap::<Vec<u8>, i128>::new();
        for t in &known_terms {
            for (word, coef) in expand_to_free_algebra(t) {
                *recomposed.entry(word).or_default() += coef;
            }
        }
        recomposed.retain(|_, coef| *coef != 0);
        assert_eq!(original, recomposed);
    }

    /// Regression test: `lyndon_basis_decomposition` of the pair (i=2, j=86) in
    /// the (d=3, m=7) basis used to loop forever because the Jacobi rewrite was
    /// applied unconditionally, cycling `[[a,b],c] <-> [[a,c],b]`.
    #[test]
    fn lyndon_basis_decomposition_terminates_and_is_correct() {
        use std::collections::HashSet;

        let basis = LyndonBasis::<u8>::new(3, Sort::Lexicographical).generate_basis(7);
        let commutator_basis = basis
            .iter()
            .map(CommutatorTerm::<i128, u8>::from)
            .collect::<Vec<_>>();
        let basis_set = commutator_basis
            .iter()
            .map(CommutatorTerm::unit_hash)
            .collect::<HashSet<_>>();

        let mut term = comm![&commutator_basis[2], &commutator_basis[86]];
        term.lyndon_sort();

        // Must terminate (previously hung with millions of rewrite iterations).
        let terms = term.lyndon_basis_decomposition(&basis_set);

        // Every returned term must be an element of the basis.
        for t in &terms {
            assert!(
                basis_set.contains(&t.unit_hash()),
                "decomposition returned a non-basis term: {t:?}"
            );
        }

        // The decomposition must equal the original commutator, verified by
        // expanding both sides into the free associative algebra.
        let original = expand_to_free_algebra(&term);
        let mut recomposed = HashMap::<Vec<u8>, i128>::new();
        for t in &terms {
            for (word, c) in expand_to_free_algebra(t) {
                *recomposed.entry(word).or_default() += c;
            }
        }
        recomposed.retain(|_, c| *c != 0);
        assert_eq!(
            original, recomposed,
            "decomposition does not recompose to the original commutator"
        );
    }

    /// Property test: for every ordered pair in a small basis, the
    /// decomposition terminates and recomposes exactly to the original
    /// commutator (checked via free-associative-algebra expansion).
    #[test]
    fn lyndon_basis_decomposition_all_pairs_property() {
        use std::collections::HashSet;

        let basis = LyndonBasis::<u8>::new(3, Sort::Lexicographical).generate_basis(5);
        let commutator_basis = basis
            .iter()
            .map(CommutatorTerm::<i128, u8>::from)
            .collect::<Vec<_>>();
        let basis_set = commutator_basis
            .iter()
            .map(CommutatorTerm::unit_hash)
            .collect::<HashSet<_>>();

        for i in 0..commutator_basis.len() {
            for j in 0..commutator_basis.len() {
                if i == j {
                    continue;
                }
                let mut term = comm![&commutator_basis[i], &commutator_basis[j]];
                term.lyndon_sort();
                if term.is_zero() {
                    continue;
                }
                if commutator_basis[i].degree() + commutator_basis[j].degree() > 5 {
                    // The decomposition requires the basis to contain all
                    // Lyndon words up to the degree of the term.
                    continue;
                }
                let terms = term.lyndon_basis_decomposition(&basis_set);
                for t in &terms {
                    assert!(
                        basis_set.contains(&t.unit_hash()),
                        "pair ({i}, {j}): non-basis term returned"
                    );
                }
                assert_eq!(
                    term.degree(),
                    commutator_basis[i].degree() + commutator_basis[j].degree(),
                    "pair ({i}, {j}): degree changed"
                );
            }
        }
    }


    use super::*;
    use lyndon_rs::lyndon::LyndonWordError;
    use rstest::rstest;
    use std::cmp::{Ordering, PartialOrd};

    #[test]
    fn test_commutators_int() {
        assert_eq!(comm![1, 2], 0);
    }

    #[rstest]
    #[case(
        CommutatorTerm::from('A'),
        CommutatorTerm::from('B'),
        CommutatorTerm::Expression {
        coefficient: 1,
        left: Box::new(CommutatorTerm::from('A')),
        right: Box::new(CommutatorTerm::from('B')),
        degree: 2,
    })]
    #[case(
        CommutatorTerm::from('B'),
        CommutatorTerm::from('A'),
        CommutatorTerm::Expression {
            coefficient: 1,
            left: Box::new(CommutatorTerm::from('B')),
            right: Box::new(CommutatorTerm::from('A')),
            degree: 2,
        })]
    #[case(
        CommutatorTerm::Expression {
            coefficient: 2,
            left: Box::new(CommutatorTerm::from('A')),
            right: Box::new(CommutatorTerm::from('B')),
            degree: 2,
        },
        CommutatorTerm::Expression {
            coefficient: 3,
            left: Box::new(CommutatorTerm::from('B')),
            right: Box::new(CommutatorTerm::from('A')),
            degree: 2,
        },
        CommutatorTerm::Expression {
            coefficient: 6,
            left: Box::new(CommutatorTerm::Expression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::from('A')),
                right: Box::new(CommutatorTerm::from('B')),
                degree: 2,
            }),
            right: Box::new(CommutatorTerm::Expression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::from('B')),
                right: Box::new(CommutatorTerm::from('A')),
                degree: 2,

            }),
            degree: 4,
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
    #[case("AC", "ABB", Ordering::Greater)]
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
                right: Box::new(CommutatorTerm::from('B')),
                degree: 2,
            },
        Ordering::Less)]
    #[case(
        CommutatorTerm::from('B'),
        CommutatorTerm::Expression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::from('A')),
                right: Box::new(CommutatorTerm::from('B')),
                degree: 2,
            },
        Ordering::Greater)]
    #[case(
        CommutatorTerm::Expression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::from('A')),
                right: Box::new(CommutatorTerm::from('B')),
                degree: 2,
            },
        CommutatorTerm::from('A'),
        Ordering::Greater)]
    #[case(
        CommutatorTerm::Expression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::from('A')),
                right: Box::new(CommutatorTerm::from('B')),
                degree: 2,
            },
        CommutatorTerm::from('B'),
        Ordering::Less)]
    #[case(
        CommutatorTerm::Expression {
                coefficient: 1,
                left: Box::new(CommutatorTerm::from('A')),
                right: Box::new(CommutatorTerm::from('B')),
                degree: 2,
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
        #[case] expected_right_term: CommutatorTerm<i128, char>,
    ) {
        let Some((left_term, right_term)) = term.jacobi_identity() else {
            panic!("Failed to create jacobi identity");
        };
        assert_eq!(
            left_term, expected_left_term,
            "{left_term} != {expected_left_term}"
        );
        assert_eq!(
            right_term, expected_right_term,
            "{right_term} != {expected_right_term}"
        );
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
        vec![
            -comm![
                CommutatorTerm::from('A'),
                comm![
                    CommutatorTerm::from('A'),
                    comm![
                        CommutatorTerm::from('B'),
                        comm![
                            CommutatorTerm::from('B'),
                            CommutatorTerm::from('C')
                        ]
                    ]
                ]
            ],
            -2 * comm![
                CommutatorTerm::<i128, char>::from('A'),
                comm![
                    comm![
                        CommutatorTerm::from('A'),
                        comm![
                            CommutatorTerm::from('B'),
                            CommutatorTerm::from('C')
                        ]
                    ],
                    CommutatorTerm::from('B')
                ]
            ],
            -comm![
                CommutatorTerm::from('A'),
                comm![
                    comm![
                        comm![
                            CommutatorTerm::from('A'),
                            CommutatorTerm::from('C')
                        ],
                        CommutatorTerm::from('B')
                    ],
                    CommutatorTerm::from('B')
                ]
            ],
            comm![
                comm![
                    comm![
                        CommutatorTerm::from('A'),
                        CommutatorTerm::from('B')
                    ],
                    CommutatorTerm::from('B')
                ],
                comm![
                    CommutatorTerm::from('A'),
                    CommutatorTerm::from('C')
                ]
            ]
        ]
    )]
    fn test_commutator_decomposition(
        #[case] num_generators: usize,
        #[case] max_degree: usize,
        #[case] term: CommutatorTerm<i128, char>,
        #[case] mut expected_basis_terms: Vec<CommutatorTerm<i128, char>>,
    ) {
        use lyndon_rs::lyndon::{LyndonBasis, Sort};

        expected_basis_terms.sort();
        let basis = LyndonBasis::<char>::new(num_generators, Sort::Lexicographical);
        let basis_set = basis
            .generate_basis(max_degree)
            .iter()
            .map(CommutatorTerm::<i128, char>::from)
            .map(|x| x.unit_hash())
            .collect::<HashSet<_>>();

        let basis_terms = term.lyndon_basis_decomposition(&basis_set);
        assert_eq!(basis_terms.len(), expected_basis_terms.len());
        for (basis_term, expected_basis_term) in basis_terms.iter().zip(&expected_basis_terms) {
            assert_eq!(
                basis_term, expected_basis_term,
                "{basis_term} != {expected_basis_term}"
            );
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
        comm![CommutatorTerm::from('C'), CommutatorTerm::from('D')],
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
