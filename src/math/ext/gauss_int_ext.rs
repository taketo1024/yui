use num_bigint::BigInt;
use num_complex::Complex;
use num_traits::{Zero, One, Signed};

use super::super::traits::*;

// TODO we should make a wrapper struct, 
// because the implementation for `inv` etc are different.

pub type GaussInt<T> = Complex<T>;

macro_rules! impl_euc_ring_gauss_int {
    ($int:ident) => {
        impl Symbol for GaussInt<$int> { 
            fn symbol() -> String {
                String::from("Z[i]")
            }
        }

        impl RingMethods for GaussInt<$int> {
            #[inline]
            fn inv(&self) -> Option<Self> {
                if self.is_unit() {
                    Some(Self::one() / self)
                } else { 
                    None
                }
            }

            #[inline]
            fn is_unit(&self) -> bool {
                (self.re.is_unit() && self.im.is_zero()) || 
                (self.re.is_zero() && self.im.is_unit())
            }

            #[inline]
            fn normalizing_unit(&self) -> Self {
                if self.re.is_positive() && !self.im.is_negative() { 
                    Self::one()
                } else if !self.re.is_positive() && self.im.is_positive() { 
                    -Self::i()
                } else if self.re.is_negative() && !self.im.is_positive() { 
                    -Self::one()
                } else if !self.re.is_negative() && self.im.is_negative() { 
                    Self::i()
                } else { 
                    Self::one()
                }
            }
        }

        decl_alg_base!(GaussInt<$int>);
        decl_add_mon!(GaussInt<$int>);
        decl_add_grp!(GaussInt<$int>);
        decl_mon!(GaussInt<$int>);
        decl_ring!(GaussInt<$int>);
        decl_euc_ring!(GaussInt<$int>);
    }
}

impl_euc_ring_gauss_int!(i32);
impl_euc_ring_gauss_int!(i64);
impl_euc_ring_gauss_int!(BigInt);

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn check_add_mon() {
        fn check<T>() where T: AddMon, for<'a> &'a T: AddMonOps<T> {}
        check::<GaussInt<i32>>();
        check::<GaussInt<i64>>();
        check::<GaussInt<BigInt>>();
    }
    
    #[test]
    fn check_add_grp() {
        fn check<T>() where T: AddGrp, for<'a> &'a T: AddGrpOps<T> {}
        check::<GaussInt<i32>>();
        check::<GaussInt<i64>>();
        check::<GaussInt<BigInt>>();
    }
    
    #[test]
    fn check_mon() {
        fn check<T>() where T: Mon, for<'a> &'a T: MonOps<T> {}
        check::<GaussInt<i32>>();
        check::<GaussInt<i64>>();
        check::<GaussInt<BigInt>>();
    }

    #[test]
    fn check_ring() {
        fn check<T>() where T: Ring, for<'a> &'a T: RingOps<T> {}
        check::<GaussInt<i32>>();
        check::<GaussInt<i64>>();
        check::<GaussInt<BigInt>>();
    }

    #[test]
    fn check_eucring() {
        fn check<T>() where T: EucRing, for<'a> &'a T: EucRingOps<T> {}
        check::<GaussInt<i32>>();
        check::<GaussInt<i64>>();
        check::<GaussInt<BigInt>>();
    }

    #[test]
    fn is_unit() { 
        assert!(GaussInt::from(1).is_unit());
        assert!(GaussInt::from(-1).is_unit());
        assert_eq!(GaussInt::from(2).is_unit(), false);
        assert!(GaussInt::new(0, 1).is_unit());
        assert!(GaussInt::new(0, -1).is_unit());
        assert_eq!(GaussInt::new(0, 2).is_unit(), false);
        assert_eq!(GaussInt::new(1, 1).is_unit(), false);
    }

    #[test]
    fn inv() { 
        // MEMO `Complex<_>::inv(&self)` is already defined.

        assert_eq!(RingMethods::inv(&GaussInt::from(1)), Some(GaussInt::from(1)));
        assert_eq!(RingMethods::inv(&GaussInt::from(-1)), Some(GaussInt::from(-1)));
        assert_eq!(RingMethods::inv(&GaussInt::from(2)), None);
        assert_eq!(RingMethods::inv(&GaussInt::new(0, 1)), Some(GaussInt::new(0, -1)));
        assert_eq!(RingMethods::inv(&GaussInt::new(0, -1)), Some(GaussInt::new(0, 1)));
        assert_eq!(RingMethods::inv(&GaussInt::new(0, 2)), None);
        assert_eq!(RingMethods::inv(&GaussInt::new(1, 1)), None);
    }

    #[test]
    fn normalizing_unit() { 
        assert_eq!(GaussInt::from(1).normalizing_unit(), GaussInt::from(1));
        assert_eq!(GaussInt::from(-1).normalizing_unit(),GaussInt::from(-1));
        assert_eq!(GaussInt::from(2).normalizing_unit(), GaussInt::from(1));
        assert_eq!(GaussInt::new(0, 1).normalizing_unit(), GaussInt::new(0, -1));
        assert_eq!(GaussInt::new(0, -1).normalizing_unit(),GaussInt::new(0, 1));
        assert_eq!(GaussInt::new(0, 2).normalizing_unit(), GaussInt::new(0, -1));
        
        assert_eq!(GaussInt::new(1, 1).normalizing_unit(), GaussInt::from(1));
        assert_eq!(GaussInt::new(-1, 1).normalizing_unit(), GaussInt::new(0, -1));
        assert_eq!(GaussInt::new(-1, -1).normalizing_unit(), GaussInt::from(-1));
        assert_eq!(GaussInt::new(1, -1).normalizing_unit(), GaussInt::new(0, 1));
    }

    #[test]
    fn divides() {
        assert!(GaussInt::from(2).divides(&GaussInt::from(4)));
        assert_eq!(GaussInt::from(3).divides(&GaussInt::from(4)), false);
        assert_eq!(GaussInt::from(0).divides(&GaussInt::from(1)), false);
        assert!(GaussInt::new(0, 2).divides(&GaussInt::from(4)));
        assert!(GaussInt::new(0, 2).divides(&GaussInt::new(0, 4)));
        assert!(GaussInt::new(1, 1).divides(&GaussInt::from(2)));
        assert_eq!(GaussInt::new(2, 1).divides(&GaussInt::from(3)), false);
    }

    #[test]
    fn gcd() {
        let (a, b) = (GaussInt::from(240), GaussInt::from(46));
        let d = GaussInt::<i32>::gcd(&a, &b);
        assert_eq!(d, GaussInt::from(2));
    }

    #[test]
    fn gcdx() {
        let (a, b) = (GaussInt::from(240), GaussInt::from(46));
        let (d, s, t) = GaussInt::<i32>::gcdx(&a, &b);
        assert_eq!(d, GaussInt::from(2));
        assert_eq!(&s * &a + &t * &b, d);
    }

    #[test]
    fn pow_mod_2() {
        assert_eq!(GaussInt::from(-1).pow_mod2(10), GaussInt::from(1));
        assert_eq!(GaussInt::from(-1).pow_mod2(-11), GaussInt::from(-1));
    }
}