use std::collections::HashMap;
use std::hash::Hash;
use std::fmt::Display;
use std::ops::{Add, Sub, Index};
use crate::math::traits::{MathElem, Ring, RingOps};
use crate::utils::format::superscript;

pub trait RModStr: Sized
where Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R> { 
    type R;

    fn zero() -> Self;
    fn rank(&self) -> usize;
    fn tors(&self) -> &Vec<Self::R>;

    fn is_free(&self) -> bool { 
        self.tors().is_empty()
    }

    fn is_zero(&self) -> bool { 
        self.rank() == 0 && self.is_free()
    }

    fn fmt_default(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use grouping_by::GroupingBy;
        let rank = self.rank();
        let tors = self.tors();

        if rank == 0 && tors.is_empty() { 
            return f.write_str("0")
        }
    
        let mut res = vec![];
        let symbol = Self::R::symbol();
    
        if rank > 1 {
            let str = format!("{}{}", symbol, superscript(rank as isize));
            res.push(str);
        } else if rank == 1 { 
            let str = format!("{}", symbol);
            res.push(str);
        }
    
        for (t, r) in tors.iter().counter(|&t| t) { 
            let str = if r > 1 { 
                format!("({}/{}){}", symbol, t, superscript(r as isize))
            } else { 
                format!("({}/{})", symbol, t)
            };
            res.push(str);
        }
    
        let str = res.join(" ⊕ ");
        f.write_str(&str) 
    }
}

#[derive(Debug, Clone)]
pub struct GenericRModStr<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    rank: usize,
    tors: Vec<R>
}

impl<R> GenericRModStr<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new(rank: usize, tors: Vec<R>) -> Self { 
        Self { rank, tors }
    }
}

impl<R> RModStr for GenericRModStr<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;

    fn zero() -> Self { 
        Self { rank: 0, tors: vec![] }
    }

    fn rank(&self) -> usize { 
        self.rank
    }

    fn tors(&self) -> &Vec<R> {
        &self.tors
    }
}

impl<R> Display for GenericRModStr<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.fmt_default(f)
    }
}

pub trait AdditiveIndex: Clone + Copy + PartialEq + Eq + Hash + Display + Add<Output = Self> + Sub<Output = Self>{}
impl <T> AdditiveIndex for T
where T: Clone + Copy + PartialEq + Eq + Hash + Display + Add<Output = Self> + Sub<Output = Self>{}

pub trait GradedRModStr: Index<Self::Index>
where 
    Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R>,
    Self::Output: RModStr<R = Self::R>
{
    type R;
    type Index: AdditiveIndex;
    type IndexRange: Iterator<Item = Self::Index>;

    fn in_range(&self, k: Self::Index) -> bool;
    fn range(&self) -> Self::IndexRange;
}

pub struct RModGrid<R, S, I, IR>
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    S: RModStr<R = R>,
    I: AdditiveIndex,
    IR: Iterator<Item = I> + Clone
{
    range: IR,
    grid: HashMap<I, S>,
    zero: S
}

impl<R, S, I, IR> RModGrid<R, S, I, IR>
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    S: RModStr<R = R>,
    I: AdditiveIndex,
    IR: Iterator<Item = I> + Clone
{
    pub fn new<F>(range: IR, f: F) -> Self
    where F: Fn(I) -> S {
        let grid = range.clone().map(|i| (i, f(i))).collect();
        let zero = S::zero();
        Self { range, grid, zero }
    }
}

impl<R, S, I, IR> Index<I> for RModGrid<R, S, I, IR>
where
    R: Ring, for<'x> &'x R: RingOps<R>,
    S: RModStr<R = R>,
    I: AdditiveIndex,
    IR: Iterator<Item = I> + Clone
{
    type Output = S;

    fn index(&self, index: I) -> &Self::Output {
        if let Some(s) = self.grid.get(&index) { 
            s
        } else { 
            &self.zero
        }
    }
}

impl<R, S, I, IR> GradedRModStr for RModGrid<R, S, I, IR>
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    S: RModStr<R = R>,
    I: AdditiveIndex,
    IR: Iterator<Item = I> + Clone
{
    type R = R;
    type Index = I;
    type IndexRange = IR;

    fn in_range(&self, k: Self::Index) -> bool {
        self.grid.contains_key(&k)
    }

    fn range(&self) -> Self::IndexRange {
        self.range.clone()
    }
}