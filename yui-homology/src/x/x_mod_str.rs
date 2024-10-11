use yui::{EucRing, EucRingOps, IndexList, Ring, RingOps};
use yui::lc::{Gen, Lc};
use yui_matrix::sparse::{SpVec, Trans};

use crate::{DisplayForGrid, RModStr, rmod_str_symbol};

// Represents a free R-module generated by a finite set of elements in X. 

#[derive(Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct XModStr<X, R>
where 
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    gens: IndexList<X>,
    rank: usize, 
    tors: Vec<R>,
    trans: Option<Trans<R>>
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

        Self { gens, rank, tors, trans }
    }

    pub fn free<Itr>(gens: Itr) -> Self 
    where Itr: IntoIterator<Item = X> {
        let gens = gens.into_iter().collect::<IndexList<X>>();
        let r = gens.len();

        Self::new(gens, r, vec![], Some(Trans::id(r))) 
    }

    pub fn zero() -> Self { 
        Self::new(IndexList::new(), 0, vec![], Some(Trans::zero()))
    }

    pub fn trans(&self) -> Option<&Trans<R>> { 
        self.trans.as_ref()
    }

    // TODO rename to `raw_gens`
    pub fn gens(&self) -> &IndexList<X> { 
        &self.gens
    }

    // MEMO Panicking is not good. 
    // Maybe this should return Option<..>, 
    // and create XFreeModStr that unwraps it. 

    pub fn gen_chain(&self, i: usize) -> Lc<X, R> { 
        if self.trans().is_none() { 
            panic!()
        }

        let n = self.dim();
        let v = SpVec::unit(n, i);

        self.as_chain(&v)
    }

    pub fn vectorize(&self, z: &Lc<X, R>) -> SpVec<R> {
        let Some(t) = self.trans() else { 
            panic!()
        };

        let n = self.gens.len();
        let v = SpVec::from_entries(n, z.iter().map(|(x, a)| { 
            let Some(i) = self.gens.index_of(x) else { 
                panic!("{x} not found in generators: {:?}", &self.gens);
            };
            (i, a.clone())
        }));

        t.forward(&v)
    }

    pub fn vectorize_euc(&self, z: &Lc<X, R>) -> SpVec<R>
    where R: EucRing, for<'x> &'x R: EucRingOps<R> {
        let r = self.rank();
        let v = self.vectorize(z);

        SpVec::from_sorted_entries(v.dim(), v.iter().map(|(i, a)| { 
            if i < r { 
                (i, a.clone())
            } else { 
                let t = &self.tors()[i - r];
                (i, a % t)
            }
        }))
    }

    pub fn as_chain(&self, v: &SpVec<R>) -> Lc<X, R> {
        let Some(t) = self.trans() else { 
            panic!()
        };

        assert_eq!(v.dim(), self.dim());

        let v = t.backward(v);

        Lc::from_iter( v.iter().map(|(i, a)| 
            (self.gens[i].clone(), a.clone())
        ) )
    }

    pub fn merge<Y>(&mut self, other: XModStr<Y, R>, reduce: bool)
    where Y: Gen { 
        self.rank = other.rank;
        self.tors = other.tors.clone();
        
        let Some(t) = self.trans.as_mut() else { 
            return;
        };

        let Some(t1) = other.trans else { 
            self.trans = None;
            return;
        };

        t.merge(t1);
        
        if reduce { 
            t.reduce()
        }
    }
}

impl<X, R> Default for XModStr<X, R>
where X: Gen, R: Ring, for<'x> &'x R: RingOps<R> {
    fn default() -> Self {
        Self::zero()
    }
}

impl<X, R> RModStr for XModStr<X, R>
where X: Gen, R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;

    fn rank(&self) -> usize {
        self.rank
    }

    fn tors(&self) -> &[R] {
        &self.tors
    }
}

impl<X, R> DisplayForGrid for XModStr<X, R>
where X: Gen, R: Ring, for<'x> &'x R: RingOps<R> {
    fn display_for_grid(&self) -> String {
        rmod_str_symbol(self.rank(), self.tors(), ".")
    }
}

#[cfg(test)]
mod tests { 
    use yui::lc::Free;

    use super::*;

    type X = Free<i32>;
    fn e(i: isize) -> X { 
        X::from(i as i32)
    }
    
    #[test]
    fn vectorize() { 
        let s = XModStr::free([e(0), e(1), e(2)]);
        
        let x = Lc::from(e(0));
        let y = Lc::from(e(1));
        let z = Lc::from(e(2));

        let v = s.vectorize(&x);
        assert_eq!(v, SpVec::unit(3, 0));

        let v = s.vectorize(&(&x - &y * 2 + &z * 3));
        assert_eq!(v, SpVec::from(vec![1,-2,3]));
    }
        
    #[test]
    fn as_chain() { 
        let s = XModStr::free([e(0), e(1), e(2)]);
        
        let x = Lc::from(e(0));
        let y = Lc::from(e(1));
        let z = Lc::from(e(2));

        let v = SpVec::from(vec![1, 2, -3]);
        let w = s.as_chain(&v);

        assert_eq!(w, &x + &y * 2 - &z * 3);
    }
}