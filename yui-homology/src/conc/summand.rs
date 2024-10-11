use yui::{EucRing, EucRingOps, IndexList, Ring, RingOps};
use yui::lc::{Gen, Lc};
use yui_matrix::sparse::{SpVec, Trans};

use crate::{DisplayForGrid, SummandTrait, rmod_str_symbol};

// Represents a free R-module generated by a finite set of elements in X. 

#[derive(Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Summand<X, R>
where 
    X: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    gens: IndexList<X>,
    rank: usize, 
    tors: Vec<R>,
    trans: Trans<R>
}

impl<X, R> Summand<X, R>
where X: Gen, R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new(gens: IndexList<X>, rank: usize, tors: Vec<R>, trans: Trans<R>) -> Self { 
        assert_eq!(trans.src_dim(), gens.len());
        assert_eq!(trans.tgt_dim(), rank + tors.len());

        Self { gens, rank, tors, trans }
    }

    pub fn free<Itr>(gens: Itr) -> Self 
    where Itr: IntoIterator<Item = X> {
        let gens = gens.into_iter().collect::<IndexList<X>>();
        let r = gens.len();

        Self::new(gens, r, vec![], Trans::id(r)) 
    }

    pub fn zero() -> Self { 
        Self::new(IndexList::new(), 0, vec![], Trans::zero())
    }

    pub fn trans(&self) -> &Trans<R> { 
        &self.trans
    }

    // TODO rename to `raw_gens`
    pub fn gens(&self) -> &IndexList<X> { 
        &self.gens
    }

    // MEMO Panicking is not good. 
    // Maybe this should return Option<..>, 
    // and create XFreeModStr that unwraps it. 

    pub fn gen_chain(&self, i: usize) -> Lc<X, R> { 
        let n = self.dim();
        let v = SpVec::unit(n, i);

        self.as_chain(&v)
    }

    pub fn vectorize(&self, z: &Lc<X, R>) -> SpVec<R> {
        let n = self.gens.len();
        let v = SpVec::from_entries(n, z.iter().map(|(x, a)| { 
            let Some(i) = self.gens.index_of(x) else { 
                panic!("{x} not found in generators: {:?}", &self.gens);
            };
            (i, a.clone())
        }));

        self.trans.forward(&v)
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
        assert_eq!(v.dim(), self.dim());

        let v = self.trans.backward(v);

        Lc::from_iter( v.iter().map(|(i, a)| 
            (self.gens[i].clone(), a.clone())
        ) )
    }

    pub fn merge<Y>(&mut self, other: Summand<Y, R>)
    where Y: Gen { 
        self.rank = other.rank;
        self.tors = other.tors.clone();
        self.trans.merge(other.trans);
        self.trans.reduce();
    }
}

impl<X, R> Default for Summand<X, R>
where X: Gen, R: Ring, for<'x> &'x R: RingOps<R> {
    fn default() -> Self {
        Self::zero()
    }
}

impl<X, R> SummandTrait for Summand<X, R>
where X: Gen, R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;

    fn rank(&self) -> usize {
        self.rank
    }

    fn tors(&self) -> &[R] {
        &self.tors
    }
}

impl<X, R> DisplayForGrid for Summand<X, R>
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
        let s = Summand::free([e(0), e(1), e(2)]);
        
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
        let s = Summand::free([e(0), e(1), e(2)]);
        
        let x = Lc::from(e(0));
        let y = Lc::from(e(1));
        let z = Lc::from(e(2));

        let v = SpVec::from(vec![1, 2, -3]);
        let w = s.as_chain(&v);

        assert_eq!(w, &x + &y * 2 - &z * 3);
    }
}