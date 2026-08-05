//! Internal storage for [`super::Lc`]: a 0/1/many specialization of a
//! `(key, coefficient)` table that avoids `HashMap` allocation for the
//! empty and single-term cases.

use std::collections::hash_map;
use rustc_hash::FxHashMap;
use crate::abst::{Ring, RingOps};

use super::lc_key::LcKey;

/// Invariants:
/// - `Single(_, r)` always has `r ≠ 0`.
/// - `Many(m)` always has `m.len() >= 2` and contains no zero values.
///
/// All mutating methods that may break these invariants either restore them
/// internally or document that [`LcData::clean`] must be called after.
#[derive(PartialEq, Eq, Clone, Debug)]
#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
pub(super) enum LcData<X, R>
where
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    Zero,
    Single(X, R),
    Many(FxHashMap<X, R>),
}

impl<X, R> Default for LcData<X, R>
where
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn default() -> Self { Self::Zero }
}

impl<X, R> LcData<X, R>
where
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    fn new_many() -> FxHashMap<X, R> {
        FxHashMap::default()
    }

    pub(super) fn len(&self) -> usize {
        match self {
            Self::Zero => 0,
            Self::Single(_, _) => 1,
            Self::Many(m) => m.len(),
        }
    }

    pub(super) fn is_empty(&self) -> bool {
        matches!(self, Self::Zero)
    }

    pub(super) fn get(&self, x: &X) -> Option<&R> {
        match self {
            Self::Zero => None,
            Self::Single(k, v) => (k == x).then_some(v),
            Self::Many(m) => m.get(x),
        }
    }

    pub(super) fn iter(&self) -> LcDataIter<'_, X, R> {
        match self {
            Self::Zero => LcDataIter::Zero,
            Self::Single(x, r) => LcDataIter::Single(Some((x, r))),
            Self::Many(m) => LcDataIter::Many(m.iter()),
        }
    }

    /// Insert/combine a `(key, coef)` pair, preserving canonical form for
    /// `Zero`/`Single` transitions. For `Many`, may leave a zero coefficient
    /// in the map; caller must invoke [`Self::clean`] afterwards.
    pub(super) fn add_pair(&mut self, x: X, r: R) {
        if r.is_zero() { return; }
        match self {
            Self::Zero => {
                *self = Self::Single(x, r);
            }
            Self::Single(x0, _) if x0 == &x => {
                let Self::Single(_, r0) = std::mem::take(self) else { unreachable!() };
                let combined = r0 + r;
                if !combined.is_zero() {
                    *self = Self::Single(x, combined);
                }
            }
            Self::Single(_, _) => {
                let Self::Single(x0, r0) = std::mem::take(self) else { unreachable!() };
                let mut m = Self::new_many();
                m.insert(x0, r0);
                m.insert(x, r);
                *self = Self::Many(m);
            }
            Self::Many(m) => {
                if let Some(v) = m.get_mut(&x) {
                    v.add_assign(r);
                } else {
                    m.insert(x, r);
                }
            }
        }
    }

    pub(super) fn add_pair_ref(&mut self, x: &X, r: &R) {
        if r.is_zero() { return; }
        // For Many: avoid cloning when the key already exists.
        if let Self::Many(m) = self {
            if let Some(v) = m.get_mut(x) {
                v.add_assign(r);
                return;
            }
        }
        self.add_pair(x.clone(), r.clone());
    }

    /// Drop zero-valued terms and downgrade `Many → Single/Zero` as needed.
    pub(super) fn clean(&mut self) {
        match self {
            Self::Zero => {}
            Self::Single(_, r) => {
                if r.is_zero() {
                    *self = Self::Zero;
                }
            }
            Self::Many(m) => {
                m.retain(|_, r| !r.is_zero());
                self.downgrade();
            }
        }
    }

    /// If `Many` has 0 or 1 entries, collapse to `Zero`/`Single`.
    fn downgrade(&mut self) {
        if let Self::Many(m) = self {
            match m.len() {
                0 => *self = Self::Zero,
                1 => {
                    let mut taken = std::mem::take(m);
                    let (x, r) = taken.drain().next().unwrap();
                    *self = Self::Single(x, r);
                }
                _ => {}
            }
        }
    }

    /// In-place coefficient mapping. Result is canonicalized.
    pub(super) fn map_coeffs_in_place<F>(&mut self, f: F)
    where F: Fn(R) -> R {
        match self {
            Self::Zero => {}
            Self::Single(_, _) => {
                let Self::Single(x, r) = std::mem::take(self) else { unreachable!() };
                let r = f(r);
                if !r.is_zero() {
                    *self = Self::Single(x, r);
                }
            }
            Self::Many(m) => {
                let taken = std::mem::take(m);
                for (k, v) in taken {
                    let v = f(v);
                    if !v.is_zero() {
                        m.insert(k, v);
                    }
                }
                self.downgrade();
            }
        }
    }
}

impl<X, R> IntoIterator for LcData<X, R>
where
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    type Item = (X, R);
    type IntoIter = LcDataIntoIter<X, R>;
    fn into_iter(self) -> Self::IntoIter {
        match self {
            Self::Zero => LcDataIntoIter::Zero,
            Self::Single(x, r) => LcDataIntoIter::Single(Some((x, r))),
            Self::Many(m) => LcDataIntoIter::Many(m.into_iter()),
        }
    }
}

/// Iterator returned by [`super::Lc::iter`].
pub enum LcDataIter<'a, X, R> {
    Zero,
    Single(Option<(&'a X, &'a R)>),
    Many(hash_map::Iter<'a, X, R>),
}

impl<'a, X, R> Iterator for LcDataIter<'a, X, R> {
    type Item = (&'a X, &'a R);
    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::Zero => None,
            Self::Single(opt) => opt.take(),
            Self::Many(it) => it.next(),
        }
    }
}

/// Owning iterator returned by `<super::Lc as IntoIterator>::into_iter`.
pub enum LcDataIntoIter<X, R> {
    Zero,
    Single(Option<(X, R)>),
    Many(hash_map::IntoIter<X, R>),
}

impl<X, R> Iterator for LcDataIntoIter<X, R> {
    type Item = (X, R);
    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::Zero => None,
            Self::Single(opt) => opt.take(),
            Self::Many(it) => it.next(),
        }
    }
}
