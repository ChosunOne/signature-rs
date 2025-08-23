use std::{
    fmt::{Debug, Display},
    hash::Hash,
    ops::{Mul, MulAssign, Neg},
};

use crate::CommutatorTerm;

/// Represents a formal indeterminate of the form `α·e₁·e₂·...·eₙ`.
///
/// A formal indeterminate is a product of symbols with a scalar coefficient.
/// This is used to represent the expansion of commutator expressions into
/// their constituent non-commutative monomial terms.
pub struct FormalIndeterminate<T, U> {
    /// The scalar coefficient multiplying the product of symbols.
    pub coefficient: U,
    /// The sequence of symbols in the product.
    pub symbols: Vec<T>,
}

impl<T, U> FormalIndeterminate<T, U> {
    /// Creates a new formal indeterminate with the given symbols and coefficient.
    #[must_use]
    pub fn new(symbols: Vec<T>, coefficient: U) -> Self {
        Self {
            coefficient,
            symbols,
        }
    }
}

impl<T: Display, U: Display> Display for FormalIndeterminate<T, U> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.coefficient)?;
        if !self.symbols.is_empty() {
            write!(f, " * ",)?;
            for symbol in &self.symbols {
                write!(f, "{symbol}")?;
            }
        }
        Ok(())
    }
}

impl<T: Debug, U: Debug> Debug for FormalIndeterminate<T, U> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("FormalIndeterminate")
            .field("coefficient", &self.coefficient)
            .field("symbols", &self.symbols)
            .finish()
    }
}

impl<T: Clone, U: Clone> Clone for FormalIndeterminate<T, U> {
    fn clone(&self) -> Self {
        Self {
            coefficient: self.coefficient.clone(),
            symbols: self.symbols.clone(),
        }
    }
}

impl<T: Hash, U: Hash> Hash for FormalIndeterminate<T, U> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.coefficient.hash(state);
        self.symbols.hash(state);
    }
}

impl<T: Clone, U: Clone + Mul<Output = U>> Mul for &FormalIndeterminate<T, U> {
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

impl<T: Clone, U: Clone + Mul<Output = U>> Mul<U> for &FormalIndeterminate<T, U> {
    type Output = FormalIndeterminate<T, U>;

    fn mul(self, rhs: U) -> Self::Output {
        FormalIndeterminate {
            coefficient: self.coefficient.clone() * rhs,
            symbols: self.symbols.clone(),
        }
    }
}

impl<T, U: Neg<Output = U>> Neg for FormalIndeterminate<T, U> {
    type Output = Self;

    fn neg(mut self) -> Self::Output {
        self.coefficient = -self.coefficient;
        self
    }
}

impl<T: Clone, U: Clone + Neg<Output = U>> Neg for &FormalIndeterminate<T, U> {
    type Output = FormalIndeterminate<T, U>;

    fn neg(self) -> Self::Output {
        let mut result = self.clone();
        result.coefficient = -result.coefficient;
        result
    }
}

impl<T: Clone, U: Clone + Mul<Output = U> + MulAssign + Neg<Output = U>> From<&CommutatorTerm<U, T>>
    for Vec<FormalIndeterminate<T, U>>
{
    fn from(value: &CommutatorTerm<U, T>) -> Self {
        match value {
            CommutatorTerm::Atom { coefficient, atom } => vec![FormalIndeterminate::<T, U> {
                coefficient: coefficient.clone(),
                symbols: vec![atom.clone()],
            }],
            CommutatorTerm::Expression {
                left: l1,
                right: r1,
                ..
            } => {
                let mut left = Self::from(&**l1);
                if let CommutatorTerm::Expression {
                    coefficient: c_l, ..
                } = &**l1
                {
                    for e_i in &mut left {
                        e_i.coefficient *= c_l.clone();
                    }
                }
                let mut right = Self::from(&**r1);
                if let CommutatorTerm::Expression {
                    coefficient: c_r, ..
                } = &**r1
                {
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

impl<T: Eq, U: Eq> PartialEq for FormalIndeterminate<T, U> {
    fn eq(&self, other: &Self) -> bool {
        self.coefficient == other.coefficient && self.symbols == other.symbols
    }
}
impl<T: Eq, U: Eq> Eq for FormalIndeterminate<T, U> {}
