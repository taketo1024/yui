//! Linear combination: a finite formal sum `Σ rᵢ · xᵢ` with `xᵢ` keys and
//! `rᵢ` coefficients in a ring `R`.
//!
//! Implements the [free `R`-module](crate::RMod) over the key set, i.e. the
//! polynomial ring viewpoint without any multiplicative structure on keys.
//!
//! Internally stored via [`super::lc_data::LcData`], which specializes the
//! empty and single-term cases (the common shape in the cobordism algebra
//! hot path of `yui-kh`) so that those paths avoid the per-entry hashmap
//! allocation. The struct [`Lc`] additionally caches an `R::zero()` so that
//! [`Lc::coeff`] can return a stable `&R` for missing keys.
//!
//! See: <https://en.wikipedia.org/wiki/Linear_combination>,
//! <https://en.wikipedia.org/wiki/Free_module>

use std::collections::HashMap;
use std::fmt::{Display, Debug};
use std::ops::{Add, AddAssign, Neg, Sub, SubAssign, Mul, MulAssign};
use itertools::Itertools;
use num_traits::Zero;
use auto_impl_ops::auto_ops;
use crate::abst::{MathType, AddMon, AddMonOps, AddGrp, AddGrpOps, Ring, RingOps, RMod, RModOps};

use super::lc_key::*;
use super::lc_data::{LcData, LcDataIter, LcDataIntoIter};

/// A linear combination `Σ rᵢ · xᵢ` with keys `X: LcKey` and coefficients in a
/// ring `R`. Stored via a private `LcData`, which specializes the empty and
/// single-term cases to avoid hashmap allocation.
#[derive(PartialEq, Eq, Clone, Default, Debug)]
#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
#[cfg_attr(feature = "serde", serde(transparent))]
pub struct Lc<X, R>
where
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    data: LcData<X, R>,
    #[cfg_attr(feature = "serde", serde(skip))]
    r_zero: R
}

impl<X, R> Lc<X, R>
where
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    pub fn new() -> Self {
        Self { data: LcData::Zero, r_zero: R::zero() }
    }

    pub fn nterms(&self) -> usize {
        self.data.len()
    }

    pub fn any_term(&self) -> Option<(&X, &R)> {
        self.iter().next()
    }

    pub fn keys(&self) -> impl Iterator<Item = &X> {
        self.iter().map(|(k, _)| k)
    }

    pub fn is_singleton(&self) -> bool {
        self.nterms() == 1 &&
        self.iter().next().unwrap().1.is_one()
    }

    pub fn as_singleton(&self) -> Option<X> {
        if !self.is_singleton() {
            None?
        }
        self.iter().next().map(|(x, _)| x.clone())
    }

    pub fn coeff(&self, x: &X) -> &R {
        self.data.get(x).unwrap_or(&self.r_zero)
    }

    pub fn iter(&self) -> LcDataIter<'_, X, R> {
        self.data.iter()
    }

    pub fn map<Y, S, F>(self, f: F) -> Lc<Y, S>
    where
        Y: LcKey,
        S: Ring, for<'x> &'x S: RingOps<S>,
        F: Fn(X, R) -> (Y, S)
    {
        self.into_iter().map(|(x, r)| f(x, r)).collect()
    }

    pub fn map_coeffs<S, F>(self, f: F) -> Lc<X, S>
    where
        S: Ring, for<'x> &'x S: RingOps<S>,
        F: Fn(R) -> S
    {
        self.map(|x, r| (x, f(r)))
    }

    pub fn map_keys<Y, F>(self, f: F) -> Lc<Y, R>
    where
        Y: LcKey,
        F: Fn(X) -> Y
    {
        self.map(|x, r| (f(x), r))
    }

    pub fn map_ref<Y, S, F>(&self, f: F) -> Lc<Y, S>
    where
        Y: LcKey,
        S: Ring, for<'x> &'x S: RingOps<S>,
        F: Fn(&X, &R) -> (Y, S)
    {
        self.iter().map(|(x, r)| f(x, r)).collect()
    }

    pub fn filter<F>(self, f: F) -> Self
    where F: Fn(&X) -> bool {
        self.into_iter().filter(|(x, _)| f(x)).collect()
    }

    pub fn filtered<F>(&self, f: F) -> Self
    where F: Fn(&X) -> bool {
        self.iter().filter_map(|(x, a)|
            if f(x) {
                Some((x.clone(), a.clone()))
            } else {
                None
            }
        ).collect()
    }

    /// Add all pairs at once. Terms may cancel along the way, so the reduced form is
    /// restored once at the end — prefer this over repeated `add_pair` in a hot loop.
    pub fn add_pairs<I>(&mut self, pairs: I)
    where I: IntoIterator<Item = (X, R)> {
        for (x, r) in pairs {
            self.data.add_pair_unreduced(x, r);
        }
        self.data.reduce();
    }

    /// Same, taking each key by reference — the coefficient is generally the cheaper
    /// of the two to clone (compare a `Cob` key in `yui-kh`), so it is passed by value.
    pub fn add_pairs_ref<'a, I>(&mut self, pairs: I)
    where I: IntoIterator<Item = (&'a X, R)>, X: 'a {
        for (x, r) in pairs {
            self.data.add_pair_ref_unreduced(x, r);
        }
        self.data.reduce();
    }

    pub fn add_pair(&mut self, rhs: (X, R)) {
        self.add_pairs([rhs]);
    }

    pub fn add_pair_ref(&mut self, rhs: (&X, R)) {
        self.add_pairs_ref([rhs]);
    }

    pub fn apply<F, Y: LcKey>(&self, f: F) -> Lc<Y, R>
    where F: Fn(&X) -> Lc<Y, R> {
        self.iter().flat_map(|(x, r)| {
            f(x).into_iter().map(move |(y, s)| {
                (y, r * &s)
            })
        }).collect()
    }

    pub fn apply_bilin<Y, Z, F>(&self, other: &Lc<Y, R>, x_map: F) -> Lc<Z, R>
    where Y: LcKey, Z: LcKey, F: Fn(&X, &Y) -> Z {
        match (&self.data, &other.data) {
            (LcData::Zero, _) | (_, LcData::Zero) => Lc::zero(),
            (LcData::Single(x, r), _) =>
                other.map_ref(|y, s| (x_map(x, y), r * s)),
            (_, LcData::Single(y, s)) =>
                self.map_ref(|x, r| (x_map(x, y), r * s)),
            (LcData::Many(_), LcData::Many(_)) => {
                let x_map = &x_map;
                let mut res = Lc::zero();
                res.add_pairs(self.iter().flat_map(|(x, r)|
                    other.iter().map(move |(y, s)| (x_map(x, y), r * s))
                ));
                res
            }
        }
    }

    pub fn sort_terms_by<F>(&self, cmp: F) -> impl Iterator<Item = (&X, &R)>
    where F: Fn(&X, &X) -> std::cmp::Ordering {
        self.iter().sorted_by(|(x, _), (y, _)| cmp(x, y))
    }

    pub fn to_string_by<F>(&self, cmp: F, descending: bool) -> String
    where F: Fn(&X, &X) -> std::cmp::Ordering {
        use crate::util::format::lc;
        if descending {
            lc( self.sort_terms_by(|x, y| cmp(x, y).reverse()) )
        } else {
            lc( self.sort_terms_by(cmp) )
        }
    }

    pub fn is_homogeneous<T, F>(&self, f: F) -> bool
    where T: PartialEq, F: Fn(&X) -> T {
        self.keys().map(f).all_equal()
    }

    pub fn homogeneous_value<T, F>(&self, f: F) -> Option<T>
    where T: PartialEq, F: Fn(&X) -> T {
        let mut iter = self.keys();
        let first = f(iter.next()?);
        if iter.all(|k| f(k) == first) { Some(first) } else { None }
    }
}

impl<X, R> From<X> for Lc<X, R>
where
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn from(x: X) -> Self {
        Self::from((x, R::one()))
    }
}

impl<X, R> From<(X, R)> for Lc<X, R>
where
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn from(value: (X, R)) -> Self {
        Self::from_iter([value])
    }
}

impl<X, R> From<HashMap<X, R>> for Lc<X, R>
where
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn from(value: HashMap<X, R>) -> Self {
        Self::from_iter(value)
    }
}

impl<X, R> FromIterator<(X, R)> for Lc<X, R>
where
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn from_iter<T: IntoIterator<Item = (X, R)>>(iter: T) -> Self {
        let mut res = Self::new();
        res.add_pairs(iter);
        res
    }
}

impl<X, R> IntoIterator for Lc<X, R>
where
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    type Item = (X, R);
    type IntoIter = LcDataIntoIter<X, R>;

    fn into_iter(self) -> Self::IntoIter {
        self.data.into_iter()
    }
}

impl<X, R> Display for Lc<X, R>
where
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.to_string_by(X::cmp, false))
    }
}

impl<X, R> Zero for Lc<X, R>
where
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn zero() -> Self {
        Self::new()
    }

    fn is_zero(&self) -> bool {
        self.data.is_empty()
    }
}

impl<X, R> Neg for Lc<X, R>
where
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        self.map_coeffs(|r| -r)
    }
}

impl<X, R> Neg for &Lc<X, R>
where
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    type Output = Lc<X, R>;

    fn neg(self) -> Self::Output {
        self.map_ref(|x, r| (x.clone(), -r))
    }
}

// Neither form clones a key that is already present. The split `auto_ops` arg-sets
// generate the four `Add` variants by rhs-ownership so the two impls don't collide.
// note: the arg-set form `auto_ops(val_val, ref_val)` is an undocumented API of `auto_impl_ops`.
#[auto_ops(val_val, ref_val)]
impl<X, R> AddAssign<Lc<X, R>> for Lc<X, R>
where
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn add_assign(&mut self, rhs: Self) {
        self.add_pairs(rhs.data);
    }
}

#[auto_ops(val_ref, ref_ref)]
impl<X, R> AddAssign<&Lc<X, R>> for Lc<X, R>
where
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn add_assign(&mut self, rhs: &Self) {
        self.add_pairs_ref(rhs.data.iter().map(|(x, r)| (x, r.clone())));
    }
}

#[auto_ops(val_val, ref_val)]
impl<X, R> SubAssign<Lc<X, R>> for Lc<X, R>
where
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn sub_assign(&mut self, rhs: Self) {
        self.add_pairs(rhs.data.into_iter().map(|(x, r)| (x, -r)));
    }
}

#[auto_ops(val_ref, ref_ref)]
impl<X, R> SubAssign<&Lc<X, R>> for Lc<X, R>
where
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn sub_assign(&mut self, rhs: &Self) {
        self.add_pairs_ref(rhs.data.iter().map(|(x, r)| (x, -r)));
    }
}

#[auto_ops]
impl<X, R> MulAssign<&R> for Lc<X, R>
where
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn mul_assign(&mut self, rhs: &R) {
        if rhs.is_one() {
            return
        }

        self.data.map_coeffs_in_place(|r| r * rhs);
    }
}

#[auto_ops]
impl<X, R> Mul for &Lc<X, R>
where
    X: LcMulKey,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    type Output = Lc<X, R>;

    fn mul(self, rhs: Self) -> Self::Output {
        self.apply_bilin(rhs, |x, y| x.mul_ref(y))
    }
}

macro_rules! impl_alg_ops {
    ($trait:ident) => {
        impl<X, R> $trait<Self> for Lc<X, R>
        where X: LcKey, R: Ring, for<'x> &'x R: RingOps<R> {}

        impl<X, R> $trait<Lc<X, R>> for &Lc<X, R>
        where X: LcKey, R: Ring, for<'x> &'x R: RingOps<R> {}
    };
}

impl_alg_ops!(AddMonOps);
impl_alg_ops!(AddGrpOps);

impl<X, R> MathType for Lc<X, R>
where
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn math_symbol() -> String {
        format!("{}<{}>", R::math_symbol(), X::math_symbol())
    }
}

impl<X, R> AddMon for Lc<X, R>
where
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>
{}

impl<X, R> AddGrp for Lc<X, R>
where
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>
{}


impl<X, R> RModOps<R, Self> for Lc<X, R>
where
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>
{}

impl<X, R> RModOps<R, Lc<X, R>> for &Lc<X, R>
where
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>
{}

impl<X, R> RMod for Lc<X, R>
where
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    type R = R;
}

#[cfg(test)]
mod tests {
    use num_traits::Zero;
    use maplit::hashmap;
    use crate::abst::{MathType, AddMon};
    use crate::lc::{AsKey, Lc};

    type X = AsKey<i32>;
    fn e(i: i32) -> X {
        X::from(i)
    }

    #[test]
    fn math_symbol() {
        type L = Lc<X, i32>;
        let symbol = L::math_symbol();
        assert_eq!(symbol, "Z<Free<i32>>");
    }

    #[test]
    fn fmt() {
        type L = Lc<X, i32>;

        let z = L::from(hashmap!{ e(1) => 1 });
        assert_eq!(z.to_string(), "<1>");

        let z = L::from(hashmap!{ e(1) => -1 });
        assert_eq!(z.to_string(), "-<1>");

        let z = L::from(hashmap!{ e(1) => 2 });
        assert_eq!(z.to_string(), "2<1>");

        let z = L::from(hashmap!{ e(1) => 1, e(2) => 1 });
        assert_eq!(z.to_string(), "<1> + <2>");

        let z = L::from(hashmap!{ e(1) => -1, e(2) => -1 });
        assert_eq!(z.to_string(), "-<1> - <2>");

        let z = L::from(hashmap!{ e(1) => 2, e(2) => 3 });
        assert_eq!(z.to_string(), "2<1> + 3<2>");

        let z = L::from(hashmap!{ e(1) => -2, e(2) => -3 });
        assert_eq!(z.to_string(), "-2<1> - 3<2>");
    }

    #[test]
    fn default() {
        type L = Lc<X, i32>;
        let z = L::default();
        assert!(z.data.is_empty());
    }

    #[test]
    fn from_singleton() {
        type L = Lc<X, i32>;
        let x = e(0);
        let z = L::from(x);
        assert_eq!(z, L::from(hashmap!{ e(0) => 1 }));
    }

    #[test]
    fn from_pair() {
        type L = Lc<X, i32>;
        let x = e(0);
        let z = L::from((x, 2));
        assert_eq!(z, L::from(hashmap!{ e(0) => 2 }));
    }

    #[test]
    fn from_iter() {
        type L = Lc<X, i32>;
        let z = L::from_iter([(e(0), 1), (e(1), 0), (e(2), 2)]);

        assert!(!z.is_zero());
        assert_eq!(z.nterms(), 2);
        assert_eq!(z.coeff(&e(0)), &1);
        assert_eq!(z.coeff(&e(2)), &2);
    }

    #[test]
    fn into_singleton() {
        type L = Lc<X, i32>;
        let z = L::from(e(0));

        assert!(z.is_singleton());
        assert_eq!(z.as_singleton(), Some(e(0)));

        let z = L::from((e(0), 2));
        assert!(!z.is_singleton());
        assert_eq!(z.as_singleton(), None);

        let z = L::from_iter([(e(0), 1), (e(1), 1)]);
        assert!(!z.is_singleton());
        assert_eq!(z.as_singleton(), None);
    }

    #[test]
    fn eq() {
        type L = Lc<X, i32>;
        let z1 = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let z2 = L::from(hashmap!{ e(2) => 2, e(1) => 1 });
        let z3 = L::from(hashmap!{ e(1) => 1 });

        assert_eq!(z1, z2);
        assert_ne!(z1, z3);
    }

    #[test]
    fn zero() {
        type L = Lc<X, i32>;
        let z = L::zero();

        assert!(z.data.is_empty());
        assert!(z.is_zero());

        let z = L::from(hashmap!{ e(1) => 1 });

        assert!(!z.data.is_empty());
        assert!(!z.is_zero());
    }

    #[test]
    fn add_pair_reduced() {
        type L = Lc<X, i32>;

        // cancelled terms must be gone the moment `add_pair` returns
        let mut z = L::from(hashmap!{ e(1) => 1, e(2) => 2, e(3) => 1 });
        z.add_pair((e(1), -1));
        assert_eq!(z, L::from(hashmap!{ e(2) => 2, e(3) => 1 }));
        assert_eq!(z.nterms(), 2);

        z.add_pair((e(2), -1));
        z.add_pair((e(3), -1));
        assert_eq!(z, L::from(hashmap!{ e(2) => 1 }));
        assert_eq!(z.nterms(), 1);
    }

    #[test]
    fn add_pairs_reduced() {
        type L = Lc<X, i32>;

        let mut z = L::from(hashmap!{ e(1) => 1, e(2) => 2, e(3) => 1 });
        z.add_pairs([(e(1), -1), (e(2), -1), (e(3), -1)]);

        assert_eq!(z, L::from(hashmap!{ e(2) => 1 }));
        assert_eq!(z.nterms(), 1);

        z.add_pairs([(e(2), -1)]);
        assert!(z.is_zero());
        assert_eq!(z.nterms(), 0);
    }

    #[test]
    fn add() {
        type L = Lc<X, i32>;
        let z1 = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let z2 = L::from(hashmap!{ e(2) => 20, e(3) => 30 });
        let w = z1 + z2;

        assert_eq!(w, L::from(hashmap!{ e(1) => 1, e(2) => 22, e(3) => 30 }));
    }

    #[test]
    fn add_ref() {
        type L = Lc<X, i32>;
        let z1 = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let z2 = L::from(hashmap!{ e(2) => 20, e(3) => 30 });
        let w = &z1 + &z2;

        assert_eq!(w, L::from(hashmap!{ e(1) => 1, e(2) => 22, e(3) => 30 }));
    }

    #[test]
    fn add_assign() {
        type L = Lc<X, i32>;
        let mut z1 = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let z2 = L::from(hashmap!{ e(2) => 20, e(3) => 30 });
        z1 += z2;

        assert_eq!(z1, L::from(hashmap!{ e(1) => 1, e(2) => 22, e(3) => 30 }));
    }

    #[test]
    fn add_assign_ref() {
        type L = Lc<X, i32>;
        let mut z1 = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let z2 = L::from(hashmap!{ e(2) => 20, e(3) => 30 });
        z1 += &z2;

        assert_eq!(z1, L::from(hashmap!{ e(1) => 1, e(2) => 22, e(3) => 30 }));
    }

    #[test]
    fn sum() {
        type L = Lc<X, i32>;
        let z1 = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let z2 = L::from(hashmap!{ e(2) => 20, e(3) => 30 });
        let z3 = L::from(hashmap!{ e(3) => 300, e(4) => 400 });
        let w  = L::sum([z1, z2, z3]);

        assert_eq!(w, L::from(hashmap!{ e(1) => 1, e(2) => 22, e(3) => 330, e(4) => 400 }));
    }

    #[test]
    fn sum_ref() {
        type L = Lc<X, i32>;
        let z1 = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let z2 = L::from(hashmap!{ e(2) => 20, e(3) => 30 });
        let z3 = L::from(hashmap!{ e(3) => 300, e(4) => 400 });
        let w  = L::sum([&z1, &z2, &z3]);

        assert_eq!(w, L::from(hashmap!{ e(1) => 1, e(2) => 22, e(3) => 330, e(4) => 400 }));
    }

    #[test]
    fn neg() {
        type L = Lc<X, i32>;
        let z = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        assert_eq!(-z, L::from(hashmap!{ e(1) => -1, e(2) => -2 }));
    }

    #[test]
    fn neg_ref() {
        type L = Lc<X, i32>;
        let z = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        assert_eq!(-(&z), L::from(hashmap!{ e(1) => -1, e(2) => -2 }));
    }

    #[test]
    fn sub() {
        type L = Lc<X, i32>;
        let z1 = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let z2 = L::from(hashmap!{ e(2) => 20, e(3) => 30 });
        let w = z1 - z2;

        assert_eq!(w, L::from(hashmap!{ e(1) => 1, e(2) => -18, e(3) => -30 }));
    }

    #[test]
    fn sub_ref() {
        type L = Lc<X, i32>;
        let z1 = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let z2 = L::from(hashmap!{ e(2) => 20, e(3) => 30 });
        let w = &z1 - &z2;

        assert_eq!(w, L::from(hashmap!{ e(1) => 1, e(2) => -18, e(3) => -30 }));
    }

    #[test]
    fn sub_assign() {
        type L = Lc<X, i32>;
        let mut z1 = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let z2 = L::from(hashmap!{ e(2) => 20, e(3) => 30 });
        z1 -= z2;

        assert_eq!(z1, L::from(hashmap!{ e(1) => 1, e(2) => -18, e(3) => -30 }));
    }

    #[test]
    fn sub_assign_ref() {
        type L = Lc<X, i32>;
        let mut z1 = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let z2 = L::from(hashmap!{ e(2) => 20, e(3) => 30 });
        z1 -= &z2;

        assert_eq!(z1, L::from(hashmap!{ e(1) => 1, e(2) => -18, e(3) => -30 }));
    }

    // The owned/borrowed `+=`/`-=` split is wired through undocumented `auto_ops` args, so
    // cross-check every generated operator form (val/ref × val/ref) and both assign forms against
    // a HashMap ground truth over Zero/Single/Many cases incl. full cancellation.
    #[test]
    fn op_forms_consistent() {
        use std::collections::HashMap;
        type L = Lc<X, i32>;

        let lc = |pairs: &[(i32, i32)]| -> L {
            L::from_iter(pairs.iter().map(|&(k, c)| (e(k), c)))
        };
        let reference = |a: &[(i32, i32)], b: &[(i32, i32)], sign: i32| -> L {
            let mut m: HashMap<i32, i32> = HashMap::new();
            for &(k, c) in a { *m.entry(k).or_default() += c; }
            for &(k, c) in b { *m.entry(k).or_default() += sign * c; }
            L::from_iter(m.into_iter().filter(|&(_, c)| c != 0).map(|(k, c)| (e(k), c)))
        };

        let cases: &[&[(i32, i32)]] = &[
            &[],
            &[(1, 5)],
            &[(1, -5)],
            &[(1, 1), (2, 2)],
            &[(2, 20), (3, 30)],
            &[(1, 3), (2, -2), (3, 7)],
            &[(1, -3), (2, 2), (3, -7)],   // negation of the previous → cancels to zero on add
            &[(1, 1), (2, 1), (3, 1), (4, 1), (5, 1)],
        ];

        for a in cases {
            for b in cases {
                let (la, lb) = (lc(a), lc(b));
                let exp_add = reference(a, b, 1);
                let exp_sub = reference(a, b, -1);

                assert_eq!(la.clone() + lb.clone(), exp_add, "Add val_val {a:?} {b:?}");
                assert_eq!(la.clone() + &lb,        exp_add, "Add val_ref {a:?} {b:?}");
                assert_eq!(&la + lb.clone(),        exp_add, "Add ref_val {a:?} {b:?}");
                assert_eq!(&la + &lb,               exp_add, "Add ref_ref {a:?} {b:?}");
                { let mut t = la.clone(); t += lb.clone(); assert_eq!(t, exp_add, "+= val {a:?} {b:?}"); }
                { let mut t = la.clone(); t += &lb;        assert_eq!(t, exp_add, "+= ref {a:?} {b:?}"); }

                assert_eq!(la.clone() - lb.clone(), exp_sub, "Sub val_val {a:?} {b:?}");
                assert_eq!(la.clone() - &lb,        exp_sub, "Sub val_ref {a:?} {b:?}");
                assert_eq!(&la - lb.clone(),        exp_sub, "Sub ref_val {a:?} {b:?}");
                assert_eq!(&la - &lb,               exp_sub, "Sub ref_ref {a:?} {b:?}");
                { let mut t = la.clone(); t -= lb.clone(); assert_eq!(t, exp_sub, "-= val {a:?} {b:?}"); }
                { let mut t = la.clone(); t -= &lb;        assert_eq!(t, exp_sub, "-= ref {a:?} {b:?}"); }

                // Borrowed operands must be untouched by the ref-rhs / ref-lhs forms.
                let _ = &la + &lb;
                let _ = &la - &lb;
                assert_eq!(la, lc(a), "lhs mutated {a:?}");
                assert_eq!(lb, lc(b), "rhs mutated {b:?}");
            }
        }
    }

    #[test]
    fn mul() {
        type L = Lc<X, i32>;
        let z = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let r = 2;
        let w = z * r;

        assert_eq!(w, L::from(hashmap!{ e(1) => 2, e(2) => 4 }));
    }

    #[test]
    fn mul_ref() {
        type L = Lc<X, i32>;
        let z = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let r = 2;
        let w = z * r;

        assert_eq!(w, L::from(hashmap!{ e(1) => 2, e(2) => 4 }));
    }

    #[test]
    fn mul_assign() {
        type L = Lc<X, i32>;
        let mut z = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let r = 2;
        z *= r;

        assert_eq!(z, L::from(hashmap!{ e(1) => 2, e(2) => 4 }));
    }

    #[test]
    fn mul_assign_ref() {
        type L = Lc<X, i32>;
        let mut z = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let r = 2;
        z *= &r;

        assert_eq!(z, L::from(hashmap!{ e(1) => 2, e(2) => 4 }));
    }

    #[test]
    fn map_coeffs() {
        type L = Lc<X, i32>;
        let z = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let w = z.map_coeffs(|a| a * 10);

        assert_eq!(w, L::from(hashmap!{ e(1) => 10, e(2) => 20 }));
    }

    #[test]
    fn map_keys() {
        type L = Lc<X, i32>;
        let z = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let w = z.map_keys(|x| e(x.0 * 10));

        assert_eq!(w, L::from(hashmap!{ e(10) => 1, e(20) => 2 }));
    }

    #[test]
    fn filter_keys() {
        type L = Lc<X, i32>;
        let z = L::from_iter( (1..10).map(|i| (e(i), i * 10)) );
        let w = z.filtered(|x| x.0 % 3 == 0 );
        assert_eq!(w, L::from(hashmap!{ e(3) => 30, e(6) => 60, e(9) => 90}))
    }

    #[test]
    #[cfg(feature = "serde")]
    fn serialize() {
        type L = Lc<X, i32>;
        let z = L::from(hashmap!{ e(1) => 1, e(2) => 2 });
        let ser = serde_json::to_string(&z).unwrap();
        let deser = serde_json::from_str::<L>(&ser).unwrap();
        assert_eq!(z, deser);
    }
}