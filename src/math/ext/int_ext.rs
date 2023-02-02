use num_bigint::BigInt;
use num_traits::Signed;
use super::super::traits::*;

pub trait IntOps<T>: EucRingOps<T> {}

pub trait Integer: EucRing + IntOps<Self> + From<i32> + Signed + num_integer::Integer
where for<'a> &'a Self: EucRingOps<Self> {}

impl<T> RingMethods for T
where T: Integer, for<'x> &'x T: IntOps<T> {
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

impl<T> DivRound for T
where T: Integer, for<'x> &'x T: IntOps<T> {
    fn div_round(&self, q: &Self) -> Self {
        use num_rational::Ratio;
        let r = Ratio::new(self.clone(), q.clone());
        r.round().to_integer()
    }
}

macro_rules! impl_ops {
    ($trait:ident, $type:ty) => {
        impl $trait<$type> for $type {}
        impl<'a> $trait<$type> for &'a $type {}
    };
}

macro_rules! impl_trait {
    ($trait:ty, $type:ty) => {
        impl $trait for $type {}
    };
    ($trait:ty, $type:ty, { $($item:item)* }) => {
        impl $trait for $type {
            $( $item )*
        }
    };
}

macro_rules! decl_integer_all {
    ($type:ident) => {
        impl_trait!(Symbol, $type, { 
            fn symbol() -> String {
                String::from("Z")
            }
        });
        impl_trait!(AlgBase, $type);

        impl_ops!(AddMonOps, $type);
        impl_trait!(AddMon, $type);

        impl_ops!(AddGrpOps, $type);
        impl_trait!(AddGrp, $type);

        impl_ops!(MonOps, $type);
        impl_trait!(Mon, $type);

        impl_ops!(RingOps, $type);
        impl_trait!(Ring, $type, {});

        impl_ops!(EucRingOps, $type);
        impl_trait!(EucRing, $type);

        impl_ops!(IntOps, $type);
        impl_trait!(Integer, $type);
    }
}

decl_integer_all!(i32);
decl_integer_all!(i64);
decl_integer_all!(BigInt);

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn check_type() {
        fn check<T>() where T: Integer, for<'a> &'a T: IntOps<T> {}
        check::<i32>();
        check::<i64>();
        check::<BigInt>();
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
}