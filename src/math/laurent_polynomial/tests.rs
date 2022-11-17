use num_traits::{One, Zero};

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
