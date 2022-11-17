use num_traits::{One, Zero};
use std::ops::{Add, Mul, Neg, Sub};
use std::{fmt, cmp};
use polynomial::Polynomial;

#[derive(Eq, PartialEq, Clone, Debug)]
struct LaurentPolynomial<T> {
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
mod tests;