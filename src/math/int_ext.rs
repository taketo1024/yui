use super::traits::*;
use super::sign::Sign;

macro_rules! impl_euc_ring_integer {
    ($type:ident) => {
        impl Symbol for $type { 
            fn symbol() -> String {
                String::from("Z")
            }
        }
        impl AlgBase for $type {}

        impl AddMonOps<$type> for $type {}
        impl<'a> AddMonOps<$type> for &'a $type {}
        impl AddMon for $type {}

        impl AddGrpOps<$type> for $type {}
        impl<'a> AddGrpOps<$type> for &'a $type {}
        impl AddGrp for $type {}

        impl MonOps<$type> for $type {}
        impl<'a> MonOps<$type> for &'a $type {}
        impl Mon for $type {}

        impl RingOps<$type> for $type {}
        impl<'a> RingOps<$type> for &'a $type {}

        impl RingMethods for $type {
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
                self == &1 || self == &-1
            }

            #[inline]
            fn normalizing_unit(&self) -> Self {
                if self >= &0 { 1 } else { -1 } 
            }
        }

        impl From<Sign> for $type { 
            #[inline]
            fn from(e: Sign) -> Self {
                if e.is_positive() { 1 } else { -1 }
            }
        }

        impl Ring for $type {}

        impl EucRingOps<$type> for $type {}
        impl<'a> EucRingOps<$type> for &'a $type {}
        impl EucRing for $type {}
    }
}

impl_euc_ring_integer!(i32);
impl_euc_ring_integer!(i64);

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn check_add_mon() {
        fn check<T>() where T: AddMon, for<'a> &'a T: AddMonOps<T> {}
        check::<i32>();
        check::<i64>();
    }
    
    #[test]
    fn check_add_grp() {
        fn check<T>() where T: AddGrp, for<'a> &'a T: AddGrpOps<T> {}
        check::<i32>();
        check::<i64>();
    }
    
    #[test]
    fn check_mon() {
        fn check<T>() where T: Mon, for<'a> &'a T: MonOps<T> {}
        check::<i32>();
        check::<i64>();
    }

    #[test]
    fn check_ring() {
        fn check<T>() where T: Ring, for<'a> &'a T: RingOps<T> {}
        check::<i32>();
        check::<i64>();
    }

    #[test]
    fn check_eucring() {
        fn check<T>() where T: EucRing, for<'a> &'a T: EucRingOps<T> {}
        check::<i32>();
        check::<i64>();
    }

    #[test]
    fn int_is_unit() { 
        assert!(1.is_unit());
        assert!((-1).is_unit());
        assert_eq!(2.is_unit(), false);
    }

    #[test]
    fn int_inv() { 
        assert_eq!(1.inv(), Some(1));
        assert_eq!((-1).inv(), Some(-1));
        assert_eq!(2.inv(), None);
    }

    #[test]
    fn int_normalizing_unit() { 
        assert_eq!(1.normalizing_unit(), 1);
        assert_eq!((-1).normalizing_unit(), -1);
        assert_eq!(2.normalizing_unit(), 1);
    }

    #[test]
    fn int_divides() {
        assert!(2.divides(&4));
        assert_eq!(3.divides(&4), false);
        assert_eq!(0.divides(&1), false);
    }

    #[test]
    fn gcd_i32() {
        let (a, b) = (240, 46);
        let d = i32::gcd(&a, &b);
        assert_eq!(d, 2);

        let (a, b) = (24, 0);
        let d = i32::gcd(&a, &b);
        assert_eq!(d, 24);

        let (a, b) = (0, 0);
        let d = i32::gcd(&a, &b);
        assert_eq!(d, 0);
    }

    #[test]
    fn gcdx_i32() {
        let (a, b) = (240, 46);
        let (d, s, t) = i32::gcdx(&a, &b);
        assert_eq!(d, 2);
        assert_eq!(&s * &a + &t * &b, d);

        let (a, b) = (24, 0);
        let (d, s, t) = i32::gcdx(&a, &b);
        assert_eq!(d, 24);
        assert_eq!(s, 1);
        assert_eq!(t, 0);

        let (a, b) = (0, 0);
        let (d, s, t) = i32::gcdx(&a, &b);
        assert_eq!(d, 0);
        assert_eq!(s, 0);
        assert_eq!(t, 0);
    }

    #[test]
    fn pow_mod_2() {
        assert_eq!((-1).pow_mod2(10), 1);
        assert_eq!((-1).pow_mod2(-11), -1);
    }
}