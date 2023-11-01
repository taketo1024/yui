use std::collections::{BTreeMap, HashSet};
use std::ops::{Add, AddAssign, Neg, SubAssign, Sub};

use auto_impl_ops::auto_ops;
use itertools::Itertools;
use num_traits::Zero;

#[derive(Clone, Default, PartialEq, Eq, Debug, Hash)]
pub struct MultiDeg<I>(BTreeMap<usize, I>) // x_0^{d_0} ... x_n^{d_n} <--> [ 0 => d_0, ..., n => d_n ]
where I: Zero;

impl<I> MultiDeg<I>
where I: Zero {
    pub fn new(mut map: BTreeMap<usize, I>) -> Self { 
        map.retain(|_, d| !d.is_zero());
        Self(map)
    }

    pub fn from_vec(degree: Vec<I>) -> Self { 
        Self::from_iter(
            degree.into_iter().enumerate()
        )
    }
    
    pub fn len(&self) -> usize { 
        self.0.len()
    }

    pub fn iter(&self) -> impl Iterator<Item = (&usize, &I)> { 
        self.0.iter()
    }

    pub fn max_index(&self) -> Option<usize> { 
        self.0.keys().max().cloned()
    }

    pub fn min_index(&self) -> Option<usize> { 
        self.0.keys().min().cloned()
    }
}

impl<I> MultiDeg<I>
where I: Clone + Zero + Ord { 
    pub fn all_leq(&self, other: &Self) -> bool {
        let i0 = HashSet::<&usize>::from_iter( self.0.keys());
        let i1 = HashSet::<&usize>::from_iter(other.0.keys());

        i0.union(&i1).all(|&&i|
            self.of(i) <= other.of(i)
        )
    }

    pub fn all_geq(&self, other: &Self) -> bool {
        other.all_leq(self)
    }
}

impl<I> From<(usize, I)> for MultiDeg<I>
where I: Zero {
    fn from(value: (usize, I)) -> Self {
        MultiDeg::from_iter([value])
    }
}

impl<I> FromIterator<(usize, I)> for MultiDeg<I> 
where I: Zero {
    fn from_iter<T: IntoIterator<Item = (usize, I)>>(iter: T) -> Self {
        Self::new(iter.into_iter().collect())
    }
}

impl<I> MultiDeg<I>
where I: Zero + Clone {
    pub fn of(&self, index: usize) -> I {
        self.0.get(&index).cloned().unwrap_or(I::zero())
    }
}

impl<I> MultiDeg<I>
where for<'x> I: Zero + Add<&'x I, Output = I> {
    pub fn total(&self) -> I { 
        self.iter().map(|(_, d)| d).fold(I::zero(), |res, d| res + d)
    }
}

impl<I> Zero for MultiDeg<I>
where for<'x> I: Clone + AddAssign<&'x I> + Zero {
    fn zero() -> Self {
        Self(BTreeMap::new())
    }

    fn is_zero(&self) -> bool {
        self.0.is_empty()
    }
}

#[auto_ops]
impl<I> AddAssign<&MultiDeg<I>> for MultiDeg<I>
where for<'x> I: AddAssign<&'x I> + Zero + Clone {
    fn add_assign(&mut self, rhs: &MultiDeg<I>) {
        let data = &mut self.0;
        for (i, d) in rhs.0.iter() { 
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
        let data = &mut self.0;
        for (i, d) in rhs.0.iter() { 
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
        let list = self.0.iter().map(|(&i, d)| 
            (i, -d)
        ).collect();
        MultiDeg(list)
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
            let i0 = HashSet::<&usize>::from_iter(self.0.keys());
            let i1 = HashSet::<&usize>::from_iter(other.0.keys());
            let indices = i0.union(&i1).sorted();
            
            for &&i in indices { 
                let c = Ord::cmp(&self.of(i), &other.of(i));
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
        assert_eq!(mdeg.of(1), -2);
        assert_eq!(mdeg.of(4), 0);
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
}