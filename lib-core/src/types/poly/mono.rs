//! Abstract monomial — a single product of variable powers, e.g. `X^a Y^b`.
//!
//! See: <https://en.wikipedia.org/wiki/Monomial>

use std::ops::{Mul, Div};
use num_traits::One;
use crate::lc::LcKey;

/// Lexicographic and graded-lex orderings on monomials.
///
/// See: <https://en.wikipedia.org/wiki/Monomial_order>
pub trait MonoOrd {
    /// Pure lexicographic order.
    fn cmp_lex(&self, other: &Self) -> std::cmp::Ordering;
    /// Graded lexicographic order (compare total degree first, break ties by lex).
    fn cmp_grlex(&self, other: &Self) -> std::cmp::Ordering;
}

/// Abstract monomial, parameterized by its degree type [`Self::Deg`].
///
/// Concrete implementors include [`Var`](super::Var), [`Var2`](super::Var2),
/// [`Var3`](super::Var3), and [`MultiVar`](super::MultiVar). With `Deg = isize`
/// the implementor represents a Laurent monomial (allowing negative powers).
pub trait Mono:
    From<Self::Deg> +
    One +
    Mul<Output = Self> +
    Div<Output = Self> +
    MonoOrd +
    LcKey
{
    /// Degree type — `usize` (ordinary) or `isize` (Laurent), or a tuple/vector
    /// for multivariate monomials.
    type Deg;

    /// The degree (or multi-degree) of this monomial.
    fn deg(&self) -> Self::Deg;

    /// `true` iff this monomial is the unit `1` (all powers zero).
    fn is_unit(&self) -> bool;

    /// Multiplicative inverse, if this monomial is a unit; otherwise `None`.
    fn inv(&self) -> Option<Self>;

    /// `true` iff `self` divides `other` (componentwise on degrees).
    fn divides(&self, other: &Self) -> bool;
}