use std::collections::BTreeMap;
use std::ops::{Add, AddAssign, Neg, SubAssign, Sub, Index};
use std::hash::Hash;

use auto_impl_ops::auto_ops;
use delegate::delegate;
use derive_more::{Display, DebugCustom};
use num_traits::Zero;

use crate::lc::OrdForDisplay;

use super::MonoOrd;

#[derive(Clone, Default, PartialEq, Eq, Hash, Display, DebugCustom)]
#[display(fmt = "{:?}", data)]
#[debug  (fmt = "{:?}", data)]
pub struct MultiDeg<I> { 
    data: BTreeMap<usize, I>, // { index => degree }
    _zero: I
}

impl<I> MultiDeg<I>
where I: Zero {
    fn new_reduced(data: BTreeMap<usize, I>) -> Self { 
        Self { data, _zero: I::zero() }
    }

    fn reduce(&mut self) { 
        self.data.retain(|_, i| !i.is_zero())
    }

    pub(crate) fn empty() -> Self { 
        Self::new_reduced(BTreeMap::new())
    }
}

impl<I> MultiDeg<I> {
    delegate! { 
        to self.data { 
            #[call(len)]
            pub fn ninds(&self) -> usize;
            pub fn iter(&self) -> impl Iterator<Item = (&usize, &I)>;
        }
    }

    pub fn indices(&self) -> impl Iterator<Item = &usize> { 
        self.data.keys()
    }

    pub fn min_index(&self) -> Option<usize> { 
        self.indices().min().cloned()
    }

    pub fn max_index(&self) -> Option<usize> { 
        self.indices().max().cloned()
    }
}

impl<I> MultiDeg<I>
where I: Zero + Ord { 
    pub fn all_leq(&self, other: &Self) -> bool {
        self.iter().all(|(&i0, d0)| { 
            d0 <= &other[i0]
        }) &&
        other.iter().all(|(&i1, d1)|
            &self[i1] <= d1
        )
    }

    pub fn all_geq(&self, other: &Self) -> bool {
        other.all_leq(self)
    }
}

impl<I> MultiDeg<I>
where I: Zero + for<'x> Add<&'x I, Output = I> {
    pub fn total(&self) -> I { 
        self.iter().map(|(_, d)| d).fold(I::zero(), |res, d| res + d)
    }
}

impl<I> From<(usize, I)> for MultiDeg<I>
where I: Zero {
    fn from(value: (usize, I)) -> Self {
        MultiDeg::from_iter([value])
    }
}

impl<I, const N: usize> From<[I; N]> for MultiDeg<I> 
where I: Zero {
    fn from(degrees: [I; N]) -> Self {
        Self::from_iter(degrees.into_iter().enumerate())
    }
}

impl<I> FromIterator<(usize, I)> for MultiDeg<I> 
where I: Zero {
    fn from_iter<T: IntoIterator<Item = (usize, I)>>(iter: T) -> Self {
        let data = iter.into_iter().filter(|(_, v)| !v.is_zero()).collect();
        Self::new_reduced(data)
    }
}

impl<I> Index<usize> for MultiDeg<I> {
    type Output = I;

    fn index(&self, i: usize) -> &Self::Output {
        self.data.get(&i).unwrap_or(&self._zero)
    }
}

impl<I> Zero for MultiDeg<I>
where I: Zero + for<'x> AddAssign<&'x I> {
    fn zero() -> Self {
        Self::empty()
    }

    fn is_zero(&self) -> bool {
        self.data.is_empty()
    }
}

#[auto_ops]
impl<I> AddAssign<&MultiDeg<I>> for MultiDeg<I>
where I: Zero + for<'x> AddAssign<&'x I> {
    fn add_assign(&mut self, rhs: &MultiDeg<I>) {
        let data = &mut self.data;
        for (i, d) in rhs.iter() { 
            if !data.contains_key(i) {
                data.insert(*i, I::zero());
            }
            let d_i = data.get_mut(i).unwrap();
            d_i.add_assign(d);
        }
        self.reduce()
    }
}

#[auto_ops]
impl<I> SubAssign<&MultiDeg<I>> for MultiDeg<I>
where I: Zero + for<'x> SubAssign<&'x I> {
    fn sub_assign(&mut self, rhs: &MultiDeg<I>) {
        let data = &mut self.data;
        for (i, d) in rhs.iter() { 
            if !data.contains_key(i) {
                data.insert(*i, I::zero());
            }
            let d_i = data.get_mut(i).unwrap();
            d_i.sub_assign(d);
        }
        self.reduce()
    }
}

impl<I> Neg for &MultiDeg<I>
where I: Zero, for<'x> &'x I: Neg<Output = I> {
    type Output = MultiDeg<I>;
    fn neg(self) -> Self::Output {
        let list = self.iter().map(|(&i, d)| 
            (i, -d)
        ).collect();
        MultiDeg::new_reduced(list)
    }
}

impl<I> MonoOrd for MultiDeg<I>
where I: Zero + Ord + for<'x> Add<&'x I, Output = I> {
    // TODO this must be fixed. 
    fn cmp_lex(&self, other: &Self) -> std::cmp::Ordering {
        let mut itr0 =  self.iter();
        let mut itr1 = other.iter();

        while let (Some((i0, d0)), Some((i1, d1))) = (itr0.next(), itr1.next()) { 
            // If i0 < i1, then x_i0 > x_i1. 
            let c = usize::cmp(i0, i1).reverse().then(I::cmp(d0, d1));
            if c.is_ne() { 
                return c;
            }
        }

        std::cmp::Ordering::Equal
    }

    fn cmp_grlex(&self, other: &Self) -> std::cmp::Ordering {
        I::cmp(&self.total(), &other.total()).then_with(|| 
            Self::cmp_lex(self, other)
        )
    }
}

impl<I> OrdForDisplay for MultiDeg<I>
where I: Zero + Ord + for<'x> Add<&'x I, Output = I> {
    fn cmp_for_display(&self, other: &Self) -> std::cmp::Ordering {
        Self::cmp_grlex(self, other)
    }
}

#[cfg(test)]
mod tests {
    use std::hash::{BuildHasher, Hasher};

    use super::*;

    #[test]
    fn reduce() { 
        let data = BTreeMap::from_iter([(1, 0), (0, 1), (7, 0), (2, 3)]);
        let mut d0 = MultiDeg{ data, _zero: 0 };

        d0.reduce();

        assert_eq!(d0.data, BTreeMap::from_iter([(0, 1), (2, 3)]));
    }

    #[test]
    fn deg() {
        let mdeg = MultiDeg::from_iter([(0, 1), (1, -2), (2, 3), (3, 0)]);
        assert_eq!(mdeg[1], -2);
        assert_eq!(mdeg[4], 0);
    }

    #[test]
    fn eq() { 
        let d1 = MultiDeg::from_iter([(2, 3), (1, -2), (3, 0), (0, 1)]);
        let d2 = MultiDeg::from_iter([(0, 1), (1, -2), (2, 3), (3, 0)]);
        let d3 = MultiDeg::from_iter([(0, 1), (1, -2), (2, 3), (3, 1)]);
        assert_eq!(d1, d2);
        assert_ne!(d1, d3);
    }

    #[test]
    fn hash() { 
        let d1 = MultiDeg::from_iter([(2, 3), (1, -2), (3, 0), (0, 1)]);
        let d2 = MultiDeg::from_iter([(0, 1), (1, -2), (2, 3), (3, 0)]);
        let d3 = MultiDeg::from_iter([(0, 1), (1, -2), (2, 3), (3, 1)]);

        let state = std::collections::hash_map::RandomState::new();
        let hash = |d: &MultiDeg<_>| -> u64 { 
            let mut hasher = state.build_hasher();
            d.hash(&mut hasher);
            hasher.finish()
        };

        assert_eq!(hash(&d1), hash(&d2));
        assert_ne!(hash(&d1), hash(&d3));
    }

    #[test]
    fn total() {
        let mdeg = MultiDeg::from_iter([(0, 1), (1, -2), (2, 3), (3, 0)]);
        assert_eq!(mdeg.total(), 2);
    }

    #[test]
    fn ord() { 
        let d0 = MultiDeg::<isize>::from_iter([]);
        let d1 = MultiDeg::from_iter([(0, 1), (1, -2), (2, 3)]); // total: 2
        let d2 = MultiDeg::from_iter([(0, 1), (1,  2), (2, 3)]); // total: 6
        let d3 = MultiDeg::from_iter([(0,-1), (1,  2), (2, 3)]); // total: 4

        assert!(MultiDeg::cmp_grlex(&d1, &d1).is_eq());
        assert!(MultiDeg::cmp_grlex(&d0, &d1).is_lt());
        assert!(MultiDeg::cmp_grlex(&d1, &d2).is_lt());
        assert!(MultiDeg::cmp_grlex(&d0, &d3).is_lt());
    }

    #[test]
    fn add() { 
        let d1 = MultiDeg::from_iter([(0, 1), (1, -2), (2, 3)]);
        let d2 = MultiDeg::from_iter([(1, 3), (2, -3), (4, 5)]);
        assert_eq!(d1 + d2, MultiDeg::from_iter([(0, 1), (1, 1), (4, 5)]))
    }

    #[test]
    fn sub() { 
        let d1 = MultiDeg::from_iter([(0, 1), (1, -2), (2, 3)]);
        let d2 = MultiDeg::from_iter([(1, 3), (2, -3), (4, 5)]);
        assert_eq!(d1 - d2, MultiDeg::from_iter([(0, 1), (1, -5), (2, 6), (4, -5)]))
    }

    #[test]
    fn sub_usize() { 
        let d1 = MultiDeg::<usize>::from_iter([(0, 1), (1, 2), (2, 3)]);
        let d2 = MultiDeg::<usize>::from_iter([(1, 2), (2, 1)]);
        assert_eq!(d1 - d2, MultiDeg::from_iter([(0, 1), (2, 2)]))
    }

    #[test]
    #[should_panic]
    fn sub_usize_panic() { 
        let d1 = MultiDeg::<usize>::from_iter([(0, 1), (1, 2), (2, 3)]);
        let d2 = MultiDeg::<usize>::from_iter([(1, 3), (2, 1)]);
        let _ = d1 - d2; // panic!
    }
}