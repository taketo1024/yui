use num_traits::{One, Zero, Pow};
use std::ops::{Add, Mul, Neg, Sub};
use std::{fmt, cmp};
use polynomial::Polynomial;

#[derive(Eq, PartialEq, Clone, Debug)]
pub struct LaurentPolynomial<T> {
    polynomial: Polynomial<T>,
    base_degree: isize
}

impl<T: Zero> LaurentPolynomial<T> {
    #[inline]
    pub fn new(data: Vec<T>, base_degree: isize) -> Self {
        // TODO truncate front
        let polynomial = Polynomial::new(data);
        Self { polynomial, base_degree }
    }

    pub fn constant(a: T) -> Self { 
        LaurentPolynomial::new(vec![a], 0)
    }
}

impl<T:One + Zero> LaurentPolynomial<T> { 
    pub fn variable() -> LaurentPolynomial<T> {
        LaurentPolynomial::new(vec![T::one()], 1)
    }
}

impl<T> LaurentPolynomial<T> where
    T: Zero + One + Eq + Neg<Output = T> + Ord + fmt::Display + Clone,
{
    pub fn pretty(&self, x: &str) -> String {
        if self.base_degree == 0 {
            self.polynomial.pretty(x)
        } else {
            format!("{}^{}({})", x, self.base_degree, self.polynomial.pretty(x)) 
        }
    }
}

impl<T> Neg for LaurentPolynomial<T> where
    T: Neg + Zero + Clone,
    <T as Neg>::Output: Zero,
{
    type Output = LaurentPolynomial<<T as Neg>::Output>;

    #[inline]
    fn neg(self) -> LaurentPolynomial<<T as Neg>::Output> {
        -&self
    }
}

impl<'a, T> Neg for &'a LaurentPolynomial<T> where
    T: Neg + Zero + Clone,
    <T as Neg>::Output: Zero,
{
    type Output = LaurentPolynomial<<T as Neg>::Output>;

    #[inline]
    fn neg(self) -> LaurentPolynomial<<T as Neg>::Output> {
        LaurentPolynomial {
            polynomial: -&self.polynomial, 
            base_degree: self.base_degree
        }
    }
}

macro_rules! forward_val_val_binop {
    (impl $imp:ident, $method:ident) => {
        impl<Lhs, Rhs> $imp<LaurentPolynomial<Rhs>> for LaurentPolynomial<Lhs>
        where
            Lhs: Zero + $imp<Rhs> + Clone,
            Rhs: Zero + Clone,
            <Lhs as $imp<Rhs>>::Output: Zero,
        {
            type Output = LaurentPolynomial<<Lhs as $imp<Rhs>>::Output>;

            #[inline]
            fn $method(self, other: LaurentPolynomial<Rhs>) -> LaurentPolynomial<<Lhs as $imp<Rhs>>::Output> {
                (&self).$method(&other)
            }
        }
    };
}

macro_rules! forward_ref_val_binop {
    (impl $imp:ident, $method:ident) => {
        impl<'a, Lhs, Rhs> $imp<LaurentPolynomial<Rhs>> for &'a LaurentPolynomial<Lhs>
        where
            Lhs: Zero + $imp<Rhs> + Clone,
            Rhs: Zero + Clone,
            <Lhs as $imp<Rhs>>::Output: Zero,
        {
            type Output = LaurentPolynomial<<Lhs as $imp<Rhs>>::Output>;

            #[inline]
            fn $method(self, other: LaurentPolynomial<Rhs>) -> LaurentPolynomial<<Lhs as $imp<Rhs>>::Output> {
                self.$method(&other)
            }
        }
    };
}

macro_rules! forward_val_ref_binop {
    (impl $imp:ident, $method:ident) => {
        impl<'a, Lhs, Rhs> $imp<&'a LaurentPolynomial<Rhs>> for LaurentPolynomial<Lhs>
        where
            Lhs: Zero + $imp<Rhs> + Clone,
            Rhs: Zero + Clone,
            <Lhs as $imp<Rhs>>::Output: Zero,
        {
            type Output = LaurentPolynomial<<Lhs as $imp<Rhs>>::Output>;

            #[inline]
            fn $method(self, other: &LaurentPolynomial<Rhs>) -> LaurentPolynomial<<Lhs as $imp<Rhs>>::Output> {
                (&self).$method(other)
            }
        }
    };
}

macro_rules! forward_all_binop {
    (impl $imp:ident, $method:ident) => {
        forward_val_val_binop!(impl $imp, $method);
        forward_ref_val_binop!(impl $imp, $method);
        forward_val_ref_binop!(impl $imp, $method);
    };
}

forward_all_binop!(impl Add, add);

impl<'a, 'b, Lhs, Rhs> Add<&'b LaurentPolynomial<Rhs>> for &'a LaurentPolynomial<Lhs> where
    Lhs: Zero + Add<Rhs> + Clone,
    Rhs: Zero + Clone,
    <Lhs as Add<Rhs>>::Output: Zero,
{
    type Output = LaurentPolynomial<<Lhs as Add<Rhs>>::Output>;

    fn add(self, other: &LaurentPolynomial<Rhs>) -> LaurentPolynomial<<Lhs as Add<Rhs>>::Output> {
        let base_degree = cmp::min(self.base_degree, other.base_degree);
        let polynomial = shift(&self.polynomial, self.base_degree - base_degree) 
                        + shift(&other.polynomial, other.base_degree - base_degree);
        LaurentPolynomial { polynomial, base_degree }
    }
}

forward_all_binop!(impl Sub, sub);

impl<'a, 'b, Lhs, Rhs> Sub<&'b LaurentPolynomial<Rhs>> for &'a LaurentPolynomial<Lhs> where
    Lhs: Zero + Sub<Rhs> + Clone,
    Rhs: Zero + Clone,
    <Lhs as Sub<Rhs>>::Output: Zero,
{
    type Output = LaurentPolynomial<<Lhs as Sub<Rhs>>::Output>;

    fn sub(self, other: &LaurentPolynomial<Rhs>) -> LaurentPolynomial<<Lhs as Sub<Rhs>>::Output> {
        let base_degree = cmp::min(self.base_degree, other.base_degree);
        let polynomial = shift(&self.polynomial, self.base_degree - base_degree) 
                        - shift(&other.polynomial, other.base_degree - base_degree);
        LaurentPolynomial { polynomial, base_degree }
    }
}

forward_all_binop!(impl Mul, mul);

impl<'a, 'b, Lhs, Rhs> Mul<&'b LaurentPolynomial<Rhs>> for &'a LaurentPolynomial<Lhs> where
    Lhs: Zero + Mul<Rhs> + Clone,
    Rhs: Zero + Clone,
    <Lhs as Mul<Rhs>>::Output: Zero,
{
    type Output = LaurentPolynomial<<Lhs as Mul<Rhs>>::Output>;

    fn mul(self, other: &LaurentPolynomial<Rhs>) -> LaurentPolynomial<<Lhs as Mul<Rhs>>::Output> {
        let base_degree = self.base_degree + other.base_degree;
        let polynomial = &self.polynomial * &other.polynomial;
        LaurentPolynomial { polynomial, base_degree }
    }
}

// TODO use macro

impl <T> Pow<i32> for &LaurentPolynomial<T> 
where
    T: Zero + One + Mul<T> + Clone
{
    type Output = LaurentPolynomial<<T as Mul<T>>::Output>;

    fn pow(self, rhs: i32) -> Self::Output {
        assert!(rhs >= 0);
        (0..rhs).fold(LaurentPolynomial::one(), |p, _| {
            p * self
        })
    }
}

impl<T: Zero + Clone> Zero for LaurentPolynomial<T> {
    #[inline]
    fn zero() -> Self {
        Self { 
            polynomial: Polynomial::zero(), 
            base_degree: 0 
        }
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self.polynomial.is_zero()
    }
}

impl<T: Zero + One + Clone> One for LaurentPolynomial<T> {
    #[inline]
    fn one() -> Self {
        Self { 
            polynomial: Polynomial::one(), 
            base_degree: 0 
        }
    }
}

fn shift<T>(polynomial: &Polynomial<T>, degree: isize) -> Polynomial<T> where
    T: Zero + Clone
{ 
    debug_assert!(degree >= 0);
    let degree = degree as usize;
    let mut data = vec![T::zero(); degree];
    for a in polynomial.data() {
        data.push(a.clone())
    }
    Polynomial::new(data)
}

#[cfg(test)]
mod tests { 
    use num_traits::{One, Zero, Pow};
    use super::LaurentPolynomial;

    #[test]
    fn new() {
        let p1 = LaurentPolynomial::new(vec![1,2,3], 0);
        let p2 = LaurentPolynomial::new(vec![1,2,3], 0);
        let p3 = LaurentPolynomial::new(vec![1,0,3], 0);
        let p4 = LaurentPolynomial::new(vec![1,2,3], -1);
        assert_eq!(p1, p2);
        assert_ne!(p1, p3);
        assert_ne!(p1, p4);
    }

    #[test]
    fn pretty() { 
        let p = LaurentPolynomial::new(vec![1,2,3], 0);
        assert_eq!(p.pretty("x"), "1+2*x+3*x^2");

        let p = LaurentPolynomial::new(vec![1,2,3], 2);
        assert_eq!(p.pretty("x"), "x^2(1+2*x+3*x^2)");

        let p = LaurentPolynomial::new(vec![1,2,3], -2);
        assert_eq!(p.pretty("x"), "x^-2(1+2*x+3*x^2)");
    }

    #[test]
    fn neg() {
        let p = LaurentPolynomial::new(vec![1,2,3], 0);
        assert_eq!(-p, LaurentPolynomial::new(vec![-1,-2,-3], 0));

        let p = LaurentPolynomial::new(vec![1,2,3], 2);
        assert_eq!(-p, LaurentPolynomial::new(vec![-1,-2,-3], 2));
    }

    #[test]
    fn add() {
        let p = LaurentPolynomial::new(vec![1,2,3], 0);
        let q = LaurentPolynomial::new(vec![-1,3], 0);
        assert_eq!(p + q, LaurentPolynomial::new(vec![0,5,3], 0));

        let p = LaurentPolynomial::new(vec![1,2,3], -1);
        let q = LaurentPolynomial::new(vec![-1,3], 1);
        assert_eq!(p + q, LaurentPolynomial::new(vec![1,2,2,3], -1));
    }

    #[test]
    fn mul() {
        let p = LaurentPolynomial::new(vec![1,2], 0);
        let q = LaurentPolynomial::new(vec![3,4], 0);
        assert_eq!(p * q, LaurentPolynomial::new(vec![3, 10, 8], 0));

        let p = LaurentPolynomial::new(vec![1,2], -3);
        let q = LaurentPolynomial::new(vec![3,4], 5);
        assert_eq!(p * q, LaurentPolynomial::new(vec![3, 10, 8], 2));
    }

    #[test]
    fn sub() {
        let p = LaurentPolynomial::new(vec![1,2,3], 0);
        let q = LaurentPolynomial::new(vec![-1,3], 0);
        assert_eq!(p - q, LaurentPolynomial::new(vec![2,-1,3], 0));

        let p = LaurentPolynomial::new(vec![1,2,3], -1);
        let q = LaurentPolynomial::new(vec![-1,3], 1);
        assert_eq!(p - q, LaurentPolynomial::new(vec![1,2,4,-3], -1));
    }

    #[test]
    fn pow() { 
        let p = LaurentPolynomial::new(vec![1,1], 0);
        assert_eq!(p.pow(0), LaurentPolynomial::one());
        assert_eq!(p.pow(1), p);
        assert_eq!(p.pow(2), LaurentPolynomial::new(vec![1,2,1], 0));
    }

    #[test]
    fn zero() {
        let p: LaurentPolynomial<i32> = LaurentPolynomial::zero();
        assert!(p.is_zero());
        assert_eq!(p, LaurentPolynomial::new(vec![], 0));
    }

    #[test]
    fn one() {
        let p: LaurentPolynomial<i32> = LaurentPolynomial::one();
        assert_eq!(p, LaurentPolynomial::new(vec![1], 0));
    }

}