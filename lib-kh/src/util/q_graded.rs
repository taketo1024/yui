//! [`QGraded`]: coefficients carrying a `q`-degree, so that the `q`-degree of a chain is
//! well defined.
//!
//! The Khovanov differential preserves `q`, and `deg h = -2`, `deg t = -4`. A coefficient that
//! mixes degrees — `1 + H`, say — leaves the chain with no `q`-degree at all, which is what
//! `is_homogeneous` reports.

use itertools::Itertools;
use num_traits::Zero;
use yui_core::abst::{Ring, RingOps};
use yui_core::num::{Ratio, FF};
use yui_core::poly::{Mono, Poly, Poly2, Var2};

use super::FastPoly;

pub trait QGraded {
    /// `false` when the element mixes `q`-degrees.
    fn is_qdeg_homogeneous(&self) -> bool;

    /// The `q`-degree. Only meaningful when [`Self::is_qdeg_homogeneous`]; zero is returned for
    /// the zero element, which is homogeneous of every degree.
    fn q_deg(&self) -> isize;
}

macro_rules! impl_scalar {
    ($($t:ty),*) => {$(
        impl QGraded for $t {
            fn is_qdeg_homogeneous(&self) -> bool {
                true
            }

            fn q_deg(&self) -> isize {
                0
            }
        }
    )*};
}

impl_scalar!(i32, i64, i128, num_bigint::BigInt);

/// `deg(a/b) = deg a - deg b`, so a graded ring of fractions grades too.
impl<T> QGraded for Ratio<T>
where T: QGraded {
    fn is_qdeg_homogeneous(&self) -> bool {
        self.numer().is_qdeg_homogeneous() && self.denom().is_qdeg_homogeneous()
    }

    fn q_deg(&self) -> isize {
        self.numer().q_deg() - self.denom().q_deg()
    }
}

impl<const P: i32> QGraded for FF<P> {
    fn is_qdeg_homogeneous(&self) -> bool {
        true
    }

    fn q_deg(&self) -> isize {
        0
    }
}

// `deg H = -2`, `deg T = -4`. A variable with no assigned degree simply has no impl.
macro_rules! impl_poly {
    ($x:literal, $deg:expr) => {
        /// `FastPoly` holds a single term, so it is always homogeneous.
        impl<R> QGraded for FastPoly<$x, R>
        where R: Zero {
            fn is_qdeg_homogeneous(&self) -> bool {
                true
            }

            fn q_deg(&self) -> isize {
                if self.coeff().is_zero() {
                    return 0;
                }
                $deg * self.deg() as isize
            }
        }

        /// A sum of terms is homogeneous only when their degrees agree.
        impl<R> QGraded for Poly<$x, R>
        where R: Ring, for<'x> &'x R: RingOps<R> {
            fn is_qdeg_homogeneous(&self) -> bool {
                self.iter().map(|(v, _)| $deg * v.deg() as isize).all_equal()
            }

            fn q_deg(&self) -> isize {
                self.iter().map(|(v, _)| $deg * v.deg() as isize).min().unwrap_or(0)
            }
        }
    };
}

impl_poly!('H', -2);
impl_poly!('T', -4);

impl<R> QGraded for Poly2<'H', 'T', R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn is_qdeg_homogeneous(&self) -> bool {
        self.iter().map(|(v, _)| ht_q_deg(v)).all_equal()
    }

    fn q_deg(&self) -> isize {
        self.iter().map(|(v, _)| ht_q_deg(v)).min().unwrap_or(0)
    }
}

fn ht_q_deg(v: &Var2<'H', 'T', usize>) -> isize {
    let (i, j) = v.multi_deg();
    -2 * i as isize - 4 * j as isize
}

#[cfg(test)]
mod tests {
    use num_traits::One;

    use super::*;

    type P  = Poly<'H', i64>;
    type Q  = Poly<'T', i64>;
    type P2 = Poly2<'H', 'T', i64>;
    type FH = FastPoly<'H', i64>;
    type FT = FastPoly<'T', i64>;

    #[test]
    fn scalars_are_degree_zero() {
        assert_eq!(3i64.q_deg(), 0);
        assert!(3i64.is_qdeg_homogeneous());
        assert_eq!(FF::<3>::from(2).q_deg(), 0);
        assert!(FF::<3>::from(2).is_qdeg_homogeneous());
    }

    #[test]
    fn fast_poly_is_always_homogeneous() {
        assert_eq!(FH::new(2, 1).q_deg(), -4);
        assert_eq!(FT::new(2, 1).q_deg(), -8);
        assert_eq!(FH::one().q_deg(), 0);
        assert!(FH::new(2, 1).is_qdeg_homogeneous());

        // zero is homogeneous of every degree; we report 0
        assert_eq!(FH::new(3, 0).q_deg(), 0);
    }

    #[test]
    fn poly_homogeneity() {
        let h = P::variable();
        assert_eq!((&h * &h).q_deg(), -4);
        assert!((&h * &h).is_qdeg_homogeneous());

        // 1 + H mixes degrees 0 and -2
        let mixed = P::one() + &h;
        assert!(!mixed.is_qdeg_homogeneous());
        assert_eq!(mixed.q_deg(), -2);

        assert_eq!((Q::variable() * Q::variable()).q_deg(), -8);
    }

    #[test]
    fn poly2_homogeneity() {
        let (h, t) = (P2::variable(0), P2::variable(1));

        assert_eq!((&h * &t).q_deg(), -6);
        assert_eq!(t.q_deg(), -4);

        // H² and T both sit in degree -4, so their sum is homogeneous
        let hh_plus_t = &h * &h + &t;
        assert!(hh_plus_t.is_qdeg_homogeneous());
        assert_eq!(hh_plus_t.q_deg(), -4);

        let mixed = &h + &t;
        assert!(!mixed.is_qdeg_homogeneous());
        assert_eq!(mixed.q_deg(), -4);
    }
}
