use std::ops::{Add, AddAssign, Neg, SubAssign, Sub, Index};
use std::hash::Hash;

use ahash::{AHashMap, AHashSet};
use auto_impl_ops::auto_ops;
use delegate::delegate;
use itertools::Itertools;
use num_traits::Zero;

#[derive(Clone, Default, PartialEq, Eq, Debug)]
pub struct MultiDeg<I> { 
    data: AHashMap<usize, I>, // { index => degree }
    _zero: I
}

impl<I> MultiDeg<I>
where I: Zero {
    fn new_raw(data: AHashMap<usize, I>) -> Self { 
        Self { data, _zero: I::zero() }
    }

    pub fn reduce(&mut self) { 
        self.data.retain(|_, i| !i.is_zero())
    }
}

impl<I> MultiDeg<I>
where I: Zero + Clone {
    pub fn reduced(&self) -> Self { 
        let mut copy = self.clone();
        copy.reduce();
        copy
    }
}

impl<I> MultiDeg<I> {
    delegate! { 
        to self.data { 
            pub fn len(&self) -> usize;
            pub fn iter(&self) -> impl Iterator<Item = (&usize, &I)>;
        }
    }

    pub fn indices(&self) -> impl Iterator<Item = &usize> {
        self.data.keys()
    }

    pub fn iter_sorted(&self) -> impl Iterator<Item = (&usize, &I)> { 
        self.data.iter().sorted_by_key(|(&i, _)| i)
    }
}

impl<I> MultiDeg<I>
where I: Zero + Ord { 
    pub fn all_leq(&self, other: &Self) -> bool {
        let i0 = AHashSet::from_iter( self.indices());
        let i1 = AHashSet::from_iter(other.indices());

        i0.union(&i1).all(|&&i|
            self[i] <= other[i]
        )
    }

    pub fn all_geq(&self, other: &Self) -> bool {
        other.all_leq(self)
    }
}

impl<I> MultiDeg<I>
where for<'x> I: Zero + Add<&'x I, Output = I> {
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
        Self::new_raw(
            iter.into_iter().filter(|(_, v)| !v.is_zero()).collect()
        )
    }
}

impl<I> Hash for MultiDeg<I>
where I: Hash {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        for (i, d) in self.iter() { 
            i.hash(state);
            d.hash(state);
        }
    }
}

impl<I> Index<usize> for MultiDeg<I> {
    type Output = I;

    fn index(&self, i: usize) -> &Self::Output {
        self.data.get(&i).unwrap_or(&self._zero)
    }
}

impl<I> Zero for MultiDeg<I>
where for<'x> I: Clone + AddAssign<&'x I> + Zero {
    fn zero() -> Self {
        Self::new_raw(AHashMap::new())
    }

    fn is_zero(&self) -> bool {
        self.data.is_empty()
    }
}

#[auto_ops]
impl<I> AddAssign<&MultiDeg<I>> for MultiDeg<I>
where for<'x> I: AddAssign<&'x I> + Zero + Clone {
    fn add_assign(&mut self, rhs: &MultiDeg<I>) {
        let data = &mut self.data;
        for (i, d) in rhs.iter() { 
            if let Some(d_i) = data.get_mut(i) { 
                d_i.add_assign(d);
            } else { 
                data.insert(i.clone(), d.clone());
            }
        }
        data.retain(|_, v| !v.is_zero())
    }
}

#[auto_ops]
impl<I> SubAssign<&MultiDeg<I>> for MultiDeg<I>
where for<'x> I: SubAssign<&'x I> + Zero + Clone {
    fn sub_assign(&mut self, rhs: &MultiDeg<I>) {
        let data = &mut self.data;
        for (i, d) in rhs.iter() { 
            if !data.contains_key(i) {
                data.insert(i.clone(), I::zero());
            }
            let d_i = data.get_mut(i).unwrap();
            d_i.sub_assign(d);
        }
        data.retain(|_, v| !v.is_zero())
    }
}

impl<I> Neg for &MultiDeg<I>
where I: Zero, for<'x> &'x I: Neg<Output = I> {
    type Output = MultiDeg<I>;
    fn neg(self) -> Self::Output {
        let list = self.iter().map(|(&i, d)| 
            (i, -d)
        ).collect();
        MultiDeg::new_raw(list)
    }
}

impl<I> PartialOrd for MultiDeg<I>
where I: Clone + Zero + Ord + for<'x> Add<&'x I, Output = I> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<I> Ord for MultiDeg<I>
where I: Clone + Zero + Ord + for<'x> Add<&'x I, Output = I> {
    // Graded lexicographic order, grlex
    // see: https://en.wikipedia.org/wiki/Monomial_order#Graded_lexicographic_order

    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        use std::cmp::*;

        Ord::cmp(&self.total(), &other.total())
        .then_with(|| { 
            let i0 = AHashSet::from_iter( self.indices());
            let i1 = AHashSet::from_iter(other.indices());
            let indices = i0.union(&i1).sorted();
            
            for &&i in indices { 
                let c = Ord::cmp(&self[i], &other[i]);
                if !c.is_eq() { 
                    return c
                }
            }
    
            Ordering::Equal
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deg() {
        let mdeg = MultiDeg::from_iter([(0, 1), (1, -2), (2, 3), (3, 0)]);
        assert_eq!(mdeg[1], -2);
        assert_eq!(mdeg[4], 0);
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

        assert!(!(d1 < d1 || d1 > d1));
        assert!(d0 < d1);
        assert!(d1 < d2);
        assert!(d0 < d3);
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

    #[test]
    fn reduced() { 
        let d0 = MultiDeg::new_raw([(1, 0), (0, 1), (7, 0), (2, 3)].into());
        let d1 = MultiDeg::new_raw([(0, 1), (2, 3)].into());

        assert_ne!(d0, d1);
        assert_eq!(d0.reduced(), d1);
    }
    
    #[test]
    fn iter_sorted() { 
        let d = MultiDeg::from_iter([(1, 2), (0, 1), (7, 1), (2, 3)]);
        let sorted = d.iter_sorted().map(|(&i, &v)| (i, v)).collect_vec();
        assert_eq!(sorted, vec![(0, 1), (1, 2), (2, 3), (7, 1)]);
    }
}