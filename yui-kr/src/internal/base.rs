use std::iter::zip;
use std::ops::RangeInclusive;

use derive_more::{Display, Debug};
use itertools::{Itertools, FoldWhile};
use num_traits::One;
use yui::{Elem, Sign, PowMod2, GetSign, Ring, RingOps};
use yui::lc::{Gen, Lc, OrdForDisplay};
use yui::poly::{PolyN, MultiVar, MonoOrd};
use yui::bitseq::{BitSeq, Bit};

pub(crate) type KRMono = MultiVar<'x', usize>;
pub(crate) type KRPoly<R> = PolyN<'x', R>;

pub type KRChain<R> = Lc<KRGen, R>;
pub type KRPolyChain<R> = Lc<KRGen, KRPoly<R>>;

#[derive(PartialEq, Eq, Hash, Default, Clone, Display, Debug)]
#[display("<{}; {}-{}>", _2, _0, _1)]
#[debug("{}", self)]
pub struct KRGen(
    pub BitSeq, // h-coords
    pub BitSeq, // v-coords
    pub KRMono
);

impl Elem for KRGen {
    fn math_symbol() -> String {
        String::from("")
    }
}

impl OrdForDisplay for KRGen {
    fn cmp_for_display(&self, other: &Self) -> std::cmp::Ordering {
        // order by priority: v > h > mono
        BitSeq::cmp(&self.1, &other.1).then_with(||
            BitSeq::cmp(&self.0, &other.0)
        ).then_with(||
            KRMono::cmp_grlex(&self.2, &other.2)
        )
    }
}

impl Gen for KRGen {}

pub(crate) fn sign_between(from: BitSeq, to: BitSeq) -> Sign { 
    use Bit::*;
    
    assert_eq!(from.len(), to.len());
    debug_assert_eq!(to.weight() - from.weight(), 1);
    
    let e = zip(from.iter(), to.iter()).fold_while(0, |c, (f, t)| 
        match (f, t) {
            (Bit0, Bit0) => FoldWhile::Continue(c),
            (Bit1, Bit1) => FoldWhile::Continue(c + 1),
            (Bit0, Bit1) => FoldWhile::Done(c),
            (Bit1, Bit0) => panic!()
        }
    ).into_inner() as u32;

    (-1).pow_mod2(e).sign()
}

pub fn combine<R>(z: KRChain<R>) -> KRPolyChain<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    z.into_iter().map(|(v, r)| {
        let w = KRGen(v.0, v.1, KRMono::one());
        let p = KRPoly::from((v.2, r));
        (w, p)
    }).collect()
}

pub fn decombine<R>(z: KRPolyChain<R>) -> KRChain<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    debug_assert!(z.iter().all(|(v, _)| v.2.is_one()));

    z.into_iter().flat_map(|(v, p)| { 
        p.into_iter().map(move |(x, a)| {
            let v = KRGen(v.0, v.1, x);
            (v, a)
        })
    }).collect()
}

pub(crate) fn extend_ends_bounded(r: RangeInclusive<isize>, d: isize, bound: RangeInclusive<isize>) -> RangeInclusive<isize> { 
    let start = isize::max(r.start() - d, *bound.start());
    let end   = isize::min(r.end()   + d, *bound.end());
    start..=end
}

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn edge_sign() {
        use Sign::*;
        
        let e = sign_between(
            BitSeq::from([0,0,0]), 
            BitSeq::from([1,0,0]));
        assert_eq!(e, Pos);

        let e = sign_between(
            BitSeq::from([0,0,0]), 
            BitSeq::from([0,1,0]));
        assert_eq!(e, Pos);

        let e = sign_between(
            BitSeq::from([1,0,0]), 
            BitSeq::from([1,1,0]));
        assert_eq!(e, Neg);

        let e = sign_between(
            BitSeq::from([0,1,0]), 
            BitSeq::from([1,1,0]));
        assert_eq!(e, Pos);
    }


}