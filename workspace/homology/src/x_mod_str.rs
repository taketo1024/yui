use std::hash::Hash;

use delegate::delegate;
use yui_core::{Ring, RingOps, IndexList};
use yui_lin_comb::{Gen, LinComb};
use yui_matrix::sparse::{SpMat, SpVec, Trans};

use crate::{DisplayForGrid, RModStr, SimpleRModStr};

// Represents a free R-module generated by a finite set of elements in X. 
#[derive(Default)]
pub struct XModStr<X, R>
where 
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    gens: IndexList<X>,
    inner: SimpleRModStr<R>
}

impl<X, R> XModStr<X, R>
where X: Gen, R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new(gens: IndexList<X>, rank: usize, tors: Vec<R>, trans: Option<Trans<R>>) -> Self { 
        if let Some(t) = &trans { 
            assert_eq!(t.src_dim(), gens.len());
            assert_eq!(t.tgt_dim(), rank + tors.len());
        } else { 
            assert_eq!(gens.len(), rank + tors.len());
        }

        let inner = SimpleRModStr::new(rank, tors, trans);
        Self::from(gens, inner)
    }

    pub fn from(gens: IndexList<X>, inner: SimpleRModStr<R>) -> Self { 
        assert_eq!(gens.len(), inner.ngens());
        Self { gens, inner }
    }

    pub fn free_from<Itr>(gens: Itr) -> Self 
    where Itr: IntoIterator<Item = X> {
        let gens = gens.into_iter().collect::<IndexList<X>>();
        let r = gens.len();

        Self::new(gens, r, vec![], None) 
    }

    pub fn zero() -> Self { 
        Self::new(IndexList::new(), 0, vec![], None)
    }

    pub fn trans(&self) -> Option<&Trans<R>> { 
        self.inner.trans()
    }

    // TODO rename to `raw_gens`
    pub fn gens(&self) -> &IndexList<X> { 
        &self.gens
    }

    pub fn ngens(&self) -> usize { 
        self.inner.ngens()
    }

    pub fn gen_chain(&self, i: usize) -> LinComb<X, R> { 
        let n = self.ngens();
        let v = SpVec::unit(n, i);
        self.as_chain(&v)
    }

    pub fn vectorize(&self, z: &LinComb<X, R>) -> SpVec<R> {
        let n = self.gens.len();
        let v = SpVec::generate(n, |set| { 
            for (x, a) in z.iter() { 
                let Some(i) = self.gens.index_of(x) else { 
                    panic!("{x} not found in generators: {:?}", &self.gens);
                };
                set(i, a.clone());
            }
        });

        if let Some(t) = self.trans() {
            t.forward(&v)
        } else { 
            v
        }
    }

    pub fn as_chain(&self, v: &SpVec<R>) -> LinComb<X, R> {
        assert_eq!(v.dim(), self.ngens());

        let v = if let Some(t) = self.trans() {
            t.backward(&v)
        } else { 
            v.clone()
        };

        let elems = v.iter().map(|(i, a)| 
            (self.gens[i].clone(), a.clone())
        );

        LinComb::from_iter(elems)
    }

    // FIXME this doesn't work when trans is some. 
    pub fn make_matrix<Y, F>(&self, to: &XModStr<Y, R>, f: F) -> SpMat<R>
    where Y: Gen, F: Fn(&X) -> Vec<(Y, R)> {
        make_matrix(&self.gens, &to.gens, f)
    }
}

impl<X, R> RModStr for XModStr<X, R>
where X: Gen, R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;

    delegate! { 
        to self.inner { 
            fn rank(&self) -> usize;
            fn tors(&self) -> &Vec<Self::R>;
        }
    }
}

impl<X, R> DisplayForGrid for XModStr<X, R>
where X: Gen, R: Ring, for<'x> &'x R: RingOps<R> {
    fn display_for_grid(&self) -> String {
        self.math_symbol()
    }
}

pub fn make_matrix<X, Y, R, F>(from: &IndexList<X>, to: &IndexList<Y>, f: F) -> SpMat<R>
where 
    X: Hash + Eq, Y: Hash + Eq, 
    R: Ring, for<'x> &'x R: RingOps<R>,
    F: Fn(&X) -> Vec<(Y, R)> 
{
    let (m, n) = (to.len(), from.len());
    SpMat::generate((m, n), |set|
        for (j, x) in from.iter().enumerate() {
            let ys = f(x);
            for (y, a) in ys {
                let i = to.index_of(&y).unwrap();
                set(i, j, a);
            }
        }
    )
}

#[cfg(test)]
mod tests { 
    use yui_lin_comb::Free;

    use super::*;

    type X = Free<i32>;
    fn e(i: isize) -> X { 
        X::from(i as i32)
    }
    
    #[test]
    fn test() { 
        let s = XModStr::free_from([e(0), e(1), e(2)]);
        
        let x = LinComb::from(e(0));
        let y = LinComb::from(e(1));
        let z = LinComb::from(e(2));

        assert_eq!(s.vectorize(&x), SpVec::unit(3, 0));
        assert_eq!(s.vectorize(&(&x - &y * 2 + &z * 3)), SpVec::from(vec![1,-2,3]));
        
        let v = SpVec::from(vec![0, 2, -3]);
        assert_eq!(s.as_chain(&v), &y * 2 - &z * 3);
    }
}