use std::collections::HashMap;
use std::fmt::Display;
use std::ops::{Add, Sub};
use std::hash::Hash;
use sprs::{CsVec, CsMat, TriMat};
use crate::math::traits::{MathElem, Ring, RingOps};
use crate::utils::collections::hashmap;
use crate::utils::format::superscript;

pub trait Summand: Sized
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

pub trait AdditiveIndex: Clone + Copy + PartialEq + Eq + Hash + Display + Add<Output = Self> + Sub<Output = Self>{}
impl <T> AdditiveIndex for T
where T: Clone + Copy + PartialEq + Eq + Hash + Display + Add<Output = Self> + Sub<Output = Self>{}

pub trait Graded {
    type Index: AdditiveIndex;
    type IndexRange: Iterator<Item = Self::Index>;

    fn in_range(&self, k: Self::Index) -> bool;
    fn range(&self) -> Self::IndexRange;
}

// -- convenient methods --

pub trait FreeGenerator: Clone + PartialEq + Eq + Hash {}
impl<T> FreeGenerator for T 
where T: Clone + PartialEq + Eq + Hash {}

pub fn vectorize<R, X>(generators: &Vec<&X>, x:&X) -> CsVec<R>
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    X: FreeGenerator
{
    vectorize_comb(generators, hashmap![x => R::one()])
}

pub fn vectorize_comb<R, X>(generators: &Vec<&X>, z:HashMap<&X, R>) -> CsVec<R>
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    X: FreeGenerator
{
    let n = generators.len();

    let mut v_ind: Vec<usize> = vec![];
    let mut v_val: Vec<R> = vec![];

    for (x, a) in z { 
        let Some(i) = generators.iter().position(|&z| x == z) else { 
            continue 
        };
        v_ind.push(i);
        v_val.push(a);
    }

    CsVec::new(n, v_ind, v_val)
}

pub fn make_matrix<R, X, F>(source: &Vec<&X>, target: &Vec<&X>, f: F) -> CsMat<R>
where 
    R: Ring, for<'x> &'x R: RingOps<R>,
    X: FreeGenerator,
    F: Fn(&X) -> Vec<(X, R)> 
{
    let (m, n) = (target.len(), source.len());

    let t_ind = target.iter()
            .enumerate()
            .map(|(i, &y)| (y, i))
            .collect::<HashMap<_, _>>();

    let mut trip = TriMat::new((m, n));

    for (j, x) in source.iter().enumerate() {
        let ys = f(x);
        for (y, a) in ys {
            if !t_ind.contains_key(&y) { continue }
            let i = t_ind[&y];
            trip.add_triplet(i, j, a);
        }
    }

    trip.to_csc()
}