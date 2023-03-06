use std::collections::BTreeMap;
use std::ops::{Add, AddAssign, Neg};

use auto_impl_ops::auto_ops;
use num_traits::Zero;

use crate::{PolyDeg, impl_polydeg_unsigned, impl_polydeg_signed};

#[derive(Clone, Default, PartialEq, Eq, Debug, Hash, PartialOrd, Ord)]
pub struct MDegree<Deg>(BTreeMap<usize, Deg>); // x_0^{d_0} ... x_n^{d_n} <--> [ 0 => d_0, ..., n => d_n ]

impl<Deg> MDegree<Deg>
where Deg: Zero {
    pub fn from_vec(degree: Vec<Deg>) -> Self { 
        Self::from_iter(
            degree.into_iter().enumerate().filter(|(_, d)| 
                !d.is_zero()
            )
        )
    }

    pub fn iter(&self) -> impl Iterator<Item = (&usize, &Deg)> { 
        self.0.iter()
    }
}

impl<Deg> From<(usize, Deg)> for MDegree<Deg> {
    fn from(value: (usize, Deg)) -> Self {
        MDegree(BTreeMap::from([value]))
    }
}

impl<Deg> FromIterator<(usize, Deg)> for MDegree<Deg> {
    fn from_iter<T: IntoIterator<Item = (usize, Deg)>>(iter: T) -> Self {
        Self(iter.into_iter().collect())
    }
}

impl<Deg> Zero for MDegree<Deg>
where for<'x> Deg: Clone + AddAssign<&'x Deg> + Zero {
    fn zero() -> Self {
        Self(BTreeMap::new())
    }

    fn is_zero(&self) -> bool {
        self.0.is_empty()
    }
}

#[auto_ops]
impl<Deg> AddAssign<&MDegree<Deg>> for MDegree<Deg>
where for<'x> Deg: Clone + AddAssign<&'x Deg> + Zero {
    fn add_assign(&mut self, rhs: &MDegree<Deg>) {
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

impl<Deg> Neg for &MDegree<Deg>
where for<'x> &'x Deg: Neg<Output = Deg> {
    type Output = MDegree<Deg>;
    fn neg(self) -> Self::Output {
        let list = self.0.iter().map(|(&i, d)| 
            (i, -d)
        ).collect();
        MDegree(list)
    }
}

impl_polydeg_unsigned!(MDegree<usize>);
impl_polydeg_signed!  (MDegree<isize>);