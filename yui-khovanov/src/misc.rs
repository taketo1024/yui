use std::collections::HashMap;
use std::ops::RangeInclusive;

use yui::lc::{Gen, Lc};
use yui::{EucRing, EucRingOps, Ring, RingOps};
use yui_homology::{isize2, Grid1, SummandTrait, Summand};
use yui_matrix::sparse::SpVec;

use crate::kh::KhChainExt;

pub fn div_vec<R>(v: &SpVec<R>, c: &R) -> Option<i32>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    v.iter().filter_map(|(_, a)|
        div(a, c)
    ).min()
}

fn div<R>(a: &R, c: &R) -> Option<i32>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    if a.is_zero() { return None }

    let mut a = a.clone();
    let mut k = 0;

    while (&a % c).is_zero() { 
        a /= c;
        k += 1;
    }

    Some(k)
}

pub fn range_of<Idx, Itr>(itr: Itr) -> RangeInclusive<Idx>
where Idx: Ord + Default + Copy, Itr: IntoIterator<Item = Idx> { 
    let (min, max) = itr.into_iter().fold(None, |res, i| { 
        if let Some((min, max)) = res { 
            if i < min { 
                Some((i, max))
            } else if max < i { 
                Some((min, i))
            } else { 
                Some((min, max))
            }
        } else {
            Some((i, i))
        }
    }).unwrap_or((Idx::default(), Idx::default()));

    min ..= max
}

pub(crate) fn collect_gen_info<X, R>(grid: &Grid1<Summand<X, R>>) -> HashMap<isize2, (usize, Vec<R>, Vec<usize>)>
where X: Gen, R: Ring, for<'x> &'x R: RingOps<R>, Lc<X, R>: KhChainExt { 
    let mut table = HashMap::new();
    let init_entry = (0, vec![], vec![]);

    for (i, h) in grid.iter() { 
        let r = h.rank();
        let t = h.tors().len();

        for k in 0..r + t { 
            let z = h.gen(k);
            let q = z.q_deg();
            let e = table.entry(isize2(i, q)).or_insert_with(|| init_entry.clone());
            if k < r { 
                e.0 += 1;
            } else { 
                e.1.push(h.tors()[k - r].clone());
            }
            e.2.push(k);
        }
    }

    table
}