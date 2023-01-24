use num_bigint::BigInt;
use num_traits::{One, Signed};
use super::int_ext::{Integer, IntOps, decl_integer};
use super::super::traits::*;


impl Symbol for BigInt { 
    fn symbol() -> String {
        String::from("Z")
    }
}

impl RingMethods for BigInt {
    #[inline]
    fn inv(&self) -> Option<Self> {
        if self.is_unit() { 
            Some(self.clone())
        } else { 
            None
        }
    }

    #[inline]
    fn is_unit(&self) -> bool {
        self.is_one() || (-self).is_one()
    }

    #[inline]
    fn normalizing_unit(&self) -> Self {
        if !self.is_negative() { 
            Self::one() 
        } else { 
            -Self::one() 
        }
    }
}

decl_alg_base!(BigInt);
decl_add_mon!(BigInt);
decl_add_grp!(BigInt);
decl_mon!(BigInt);
decl_ring!(BigInt);
decl_euc_ring!(BigInt);
decl_integer!(BigInt);

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn check_add_mon() {
        fn check<T>() where T: AddMon, for<'a> &'a T: AddMonOps<T> {}
        check::<BigInt>();
    }
    
    #[test]
    fn check_add_grp() {
        fn check<T>() where T: AddGrp, for<'a> &'a T: AddGrpOps<T> {}
        check::<BigInt>();
    }
    
    #[test]
    fn check_mon() {
        fn check<T>() where T: Mon, for<'a> &'a T: MonOps<T> {}
        check::<BigInt>();
    }

    #[test]
    fn check_ring() {
        fn check<T>() where T: Ring, for<'a> &'a T: RingOps<T> {}
        check::<BigInt>();
    }

    #[test]
    fn check_eucring() {
        fn check<T>() where T: EucRing, for<'a> &'a T: EucRingOps<T> {}
        check::<BigInt>();
    }

    #[test]
    fn is_unit() { 
        assert!(BigInt::from(1).is_unit());
        assert!(BigInt::from(-1).is_unit());
        assert_eq!(BigInt::from(2).is_unit(), false);
    }

    #[test]
    fn inv() { 
        assert_eq!(BigInt::from(1).inv(), Some(BigInt::from(1)));
        assert_eq!(BigInt::from(-1).inv(), Some(BigInt::from(-1)));
        assert_eq!(BigInt::from(2).inv(), None);
    }

    #[test]
    fn normalizing_unit() { 
        assert_eq!(BigInt::from(1).normalizing_unit(), BigInt::from(1));
        assert_eq!(BigInt::from(-1).normalizing_unit(), BigInt::from(-1));
        assert_eq!(BigInt::from(2).normalizing_unit(), BigInt::from(1));
    }

    #[test]
    fn divides() {
        assert!(BigInt::from(2).divides(&BigInt::from(4)));
        assert_eq!(BigInt::from(3).divides(&BigInt::from(4)), false);
        assert_eq!(BigInt::from(0).divides(&BigInt::from(1)), false);
    }

    #[test]
    fn gcd() {
        let (a, b) = (BigInt::from(240), BigInt::from(46));
        let d = BigInt::gcd(&a, &b);
        assert_eq!(d, BigInt::from(2));

        let (a, b) = (BigInt::from(24), BigInt::from(0));
        let d = BigInt::gcd(&a, &b);
        assert_eq!(d, BigInt::from(24));

        let (a, b) = (BigInt::from(0), BigInt::from(0));
        let d = BigInt::gcd(&a, &b);
        assert_eq!(d, BigInt::from(0));
    }

    #[test]
    fn gcdx() {
        let (a, b) = (BigInt::from(240), BigInt::from(46));
        let (d, s, t) = BigInt::gcdx(&a, &b);
        assert_eq!(d, BigInt::from(2));
        assert_eq!(&s * &a + &t * &b, d);

        let (a, b) = (BigInt::from(24), BigInt::from(0));
        let (d, s, t) = BigInt::gcdx(&a, &b);
        assert_eq!(d, BigInt::from(24));
        assert_eq!(s, BigInt::from(1));
        assert_eq!(t, BigInt::from(0));

        let (a, b) = (BigInt::from(0), BigInt::from(0));
        let (d, s, t) = BigInt::gcdx(&a, &b);
        assert_eq!(d, BigInt::from(0));
        assert_eq!(s, BigInt::from(0));
        assert_eq!(t, BigInt::from(0));
    }

    #[test]
    fn pow_mod_2() {
        assert_eq!(BigInt::from(-1).pow_mod2(10), BigInt::from(1));
        assert_eq!(BigInt::from(-1).pow_mod2(-11), BigInt::from(-1));
    }
}