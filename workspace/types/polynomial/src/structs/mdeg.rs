use std::collections::BTreeMap;
use std::ops::{Add, AddAssign, Neg};

use auto_impl_ops::auto_ops;
use num_traits::Zero;

use crate::{MonoDeg, impl_deg_unsigned, impl_deg_signed};

#[derive(Clone, Default, PartialEq, Eq, Debug, Hash, PartialOrd, Ord)]
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
    pub fn deg(&self, index: usize) -> I {
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

impl_deg_unsigned!(MultiDeg<usize>);
impl_deg_signed!  (MultiDeg<isize>);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deg() {
        let mdeg = MultiDeg::from_iter([(0, 1), (1, -2), (2, 3), (3, 0)]);
        assert_eq!(mdeg.deg(1), -2);
        assert_eq!(mdeg.deg(4), 0);
    }

    #[test]
    fn total() {
        let mdeg = MultiDeg::from_iter([(0, 1), (1, -2), (2, 3), (3, 0)]);
        assert_eq!(mdeg.total(), 2);
    }
}