use std::iter::Sum;
use std::ops::Add;

use itertools::Itertools;

use yui::{Ring, RingOps};
use yui_matrix::sparse::{Trans, SpVec, SpMat};

use crate::DisplayForGrid;

pub trait RModStr
where Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R> {
    type R;

    fn rank(&self) -> usize;
    fn tors(&self) -> &[Self::R];

    fn dim(&self) -> usize { // not a good name...
        self.rank() + self.tors().len()
    }

    fn is_zero(&self) -> bool { 
        self.rank() == 0 && self.is_free()
    }

    fn is_free(&self) -> bool { 
        self.tors().is_empty()
    }

    fn math_symbol(&self) -> String { 
        rmod_str_symbol(self.rank(), self.tors(), "0")
    }
}

pub fn rmod_str_symbol<R>(rank: usize, tors: &[R], dflt: &str) -> String
where R: Ring, for<'x> &'x R: RingOps<R> {
    use yui::util::format::superscript;

    let tors = tors.iter()
        .into_group_map_by(|r| r.to_string())
        .into_iter().map(|(k, list)| (k, list.len()))
        .collect_vec();

    if rank == 0 && tors.is_empty() { 
        return dflt.to_string()
    }

    let mut res = vec![];
    let symbol = R::math_symbol();

    if rank > 1 {
        let str = format!("{}{}", symbol, superscript(rank as isize));
        res.push(str);
    } else { 
        let str = symbol.to_string();
        res.push(str);
    }
    
    for (t, r) in tors.iter() { 
        let str = if r > &1 { 
            format!("({}/{}){}", symbol, t, superscript(*r as isize))
        } else { 
            format!("({}/{})", symbol, t)
        };
        res.push(str);
    }

    
    res.join(" ⊕ ")
}

#[derive(Clone, Debug)]
pub struct SimpleRModStr<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    rank: usize, 
    tors: Vec<R>,
    trans: Option<Trans<R>>
}

impl<R> SimpleRModStr<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new(rank: usize, tors: Vec<R>, trans: Option<Trans<R>>) -> Self { 
        if let Some(t) = &trans { 
            assert_eq!(rank + tors.len(), t.tgt_dim());
        }
        Self { rank, tors, trans }
    }

    pub fn free(rank: usize) -> Self { 
        Self::new(rank, vec![], Some(Trans::id(rank)))
    }

    pub fn zero() -> Self { 
        Self::new(0, vec![], Some(Trans::zero()))
    }

    pub fn trans(&self) -> Option<&Trans<R>> { 
        self.trans.as_ref()
    }

    // The vector representing the i-th generator.

    pub fn gen_vec(&self, i: usize) -> SpVec<R> {
        let Some(t) = &self.trans else { 
            panic!()
        };

        assert!(i < self.dim());

        t.backward_mat().col_vec(i)
    }

    pub fn compose(&self, other: &SimpleRModStr<R>) -> SimpleRModStr<R> { 
        let rank = other.rank;
        let tors = other.tors.clone();

        if let Some(t0) = &self.trans { 
            if let Some(t1) = &other.trans { 
                let t = t0.compose(t1);
                return Self::new(rank, tors, Some(t))
            }
        }
        Self::new(rank, tors, None)
    }
}

impl<R> Default for SimpleRModStr<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn default() -> Self {
        Self::zero()
    }
}

impl<R> RModStr for SimpleRModStr<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;

    fn rank(&self) -> usize {
        self.rank
    }

    fn tors(&self) -> &[R] {
        &self.tors
    }
}

impl<R> DisplayForGrid for SimpleRModStr<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn display_for_grid(&self) -> String {
        rmod_str_symbol(self.rank(), self.tors(), ".")
    }
}

// direct sum
impl<'a, R> Sum<&'a SimpleRModStr<R>> for SimpleRModStr<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn sum<I: Iterator<Item = &'a SimpleRModStr<R>>>(iter: I) -> Self {
        let mut rank = 0;
        let mut tors = vec![];
        let mut src_rank = 0;
        let mut f_entries = vec![];
        let mut b_entries = vec![];
        let mut trans_ok = true;

        iter.for_each(|s| { 
            if trans_ok {
                if let Some(t) = s.trans() { 
                    f_entries.extend( t.forward_mat() .iter().map(|(i, j, a)| (i + rank, j + src_rank, a.clone())) );
                    b_entries.extend( t.backward_mat().iter().map(|(i, j, a)| (i + src_rank, j + rank, a.clone())) );
                    src_rank += t.src_dim();
                } else { 
                    trans_ok = false;
                }
            }

            rank += s.rank;
            tors.append(&mut s.tors.clone());
        });

        let trans = if trans_ok { 
            let f = SpMat::from_entries((rank, src_rank), f_entries);
            let b = SpMat::from_entries((src_rank, rank), b_entries);
            Some(Trans::new(f, b))
        } else { 
            None
        };

        SimpleRModStr::new(rank, tors, trans)
    }
}

impl<'a, 'b, R> Add<&'b SimpleRModStr<R>> for &'a SimpleRModStr<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type Output = SimpleRModStr<R>;

    fn add(self, other: &'b SimpleRModStr<R>) -> Self::Output {
        [self, other].into_iter().sum()
    }
}