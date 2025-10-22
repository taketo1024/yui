use num_bigint::BigInt;
use num_traits::{One, Signed, ToPrimitive, FromPrimitive};
use crate::*;

pub trait IntOps<T = Self>: EucRingOps<T> {}

pub trait Integer: EucRing + IntOps + Signed + PartialOrd + Ord + FromPrimitive + ToPrimitive
where for<'a> &'a Self: EucRingOps<Self> {}

impl<T> DivRound for T
where T: Integer, for<'x> &'x T: IntOps<T> {
    fn div_round(&self, q: &Self) -> Self {
        let a = self.to_f64().unwrap();
        let b = q.to_f64().unwrap();
        let r = (a / b).round();
        Self::from_f64(r).unwrap()
    }
}

macro_rules! impl_ops {
    ($trait:ident, $type:ty) => {
        impl $trait for $type {}
        impl<'a> $trait<$type> for &'a $type {}
    };
}

macro_rules! impl_integer {
    ($type:ident) => {
        impl_ops!(AddMonOps, $type);
        impl_ops!(AddGrpOps, $type);
        impl_ops!(MonOps, $type);
        impl_ops!(RingOps, $type);
        impl_ops!(EucRingOps, $type);
        impl_ops!(IntOps, $type);

        impl Elem for $type {
            fn math_symbol() -> String { 
                String::from("Z")
            }
        }
        
        impl AddMon for $type {}
        impl AddGrp for $type {}
        impl Mon for $type {}
        impl Ring for $type {
            fn inv(&self) -> Option<Self> {
                if self.is_unit() { 
                    Some(self.clone())
                } else { 
                    None
                }
            }
        
            fn is_unit(&self) -> bool {
                self.is_one() || (-self).is_one()
            }
        
            fn normalizing_unit(&self) -> Self {
                if !self.is_negative() { 
                    Self::one() 
                } else { 
                    -Self::one() 
                }
            }

            fn c_weight(&self) -> f64 {
                self.abs().to_f64().unwrap()
            }
        }

        impl EucRing for $type {
            fn gcd(x: &Self, y: &Self) -> Self {
                num_integer::Integer::gcd(x, y)
            }

            fn gcdx(x: &Self, y: &Self) -> (Self, Self, Self) {
                let num_integer::ExtendedGcd{ gcd: d, x: s, y: t } = num_integer::Integer::extended_gcd(x, y);
                (d, s, t)
            }

            fn lcm(x: &Self, y: &Self) -> Self {
                num_integer::Integer::lcm(x, y)
            }
        }

        impl Integer for $type {}
    }
}

impl_integer!(i32);
impl_integer!(i64);
impl_integer!(i128);
impl_integer!(BigInt);

#[cfg(feature = "tex")] 
mod tex {
    use crate::tex::TeX;
    use num_bigint::BigInt;

    macro_rules! impl_tex {
        ($type:ident) => {
            impl TeX for $type {
                fn tex_math_symbol() -> String { 
                    String::from("\\mathbb{Z}")
                }
                fn tex_string(&self) -> String {
                    self.to_string()
                }
            }
        }
    }

    impl_tex!(i32);
    impl_tex!(i64);
    impl_tex!(i128);
    impl_tex!(BigInt);        
}

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn check_type() {
        fn check<T>() where T: Integer, for<'a> &'a T: IntOps<T> {}
        check::<i32>();
        check::<i64>();
        check::<i128>();
        check::<BigInt>();
    }

    #[test]
    fn int_is_unit() { 
        assert!(1.is_unit());
        assert!((-1).is_unit());
        assert!(!2.is_unit());
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
        assert!(!3.divides(&4));
        assert!(!0.divides(&1));
    }

    #[test]
    fn gcd_i32() {
        let (a, b) = (240, 46);
        let d = i32::gcd(&a, &b);
        assert_eq!(d, 2);

        let (a, b) = (24, 0);
        let d = i32::gcd(&a, &b);
        assert_eq!(d, 24);

        let (a, b) = (0, -24);
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
        assert_eq!(s * a + t * b, d);

        let (a, b) = (24, 0);
        let (d, s, t) = i32::gcdx(&a, &b);
        assert_eq!(d, 24);
        assert_eq!(s * a + t * b, d);

        let (a, b) = (0, 0);
        let (d, s, t) = i32::gcdx(&a, &b);
        assert_eq!(d, 0);
        assert_eq!(s * a + t * b, d);
    }

    #[test]
    fn div_round() { 
        assert_eq!(12.div_round(&5), 2);
        assert_eq!(13.div_round(&5), 3);
        assert_eq!((-12).div_round(&5), -2);
        assert_eq!((-13).div_round(&5), -3);
    }

    #[cfg(feature = "tex")]
    #[test]
    fn tex() { 
        use crate::tex::*;
        assert_eq!(i32::tex_math_symbol(), "\\mathbb{Z}");
        assert_eq!((-2).tex_string(), "-2");
    }
}