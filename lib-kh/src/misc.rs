use std::collections::HashMap;
use std::ops::RangeInclusive;

use yui_core::lc::{Gen, Lc};
use yui_core::{EucRing, EucRingOps, Ring, RingOps};
use yui_homology::{isize2, Grid1, Grid2, Summand, SummandTrait};
use yui_matrix::sparse::SpVec;

use crate::kh::KhChainExt;
use cartesian::cartesian;

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

pub(crate) fn make_gen_grid<X, R>(grid: &Grid1<Summand<X, R>>) -> Grid2<Summand<X, R>> 
where X: Gen, R: Ring, for<'x> &'x R: RingOps<R>, Lc<X, R>: KhChainExt { 
    let info = collect_gen_info(grid);

    let h_range = range_of(info.keys().map(|i| i.0));
    let q_range = range_of(info.keys().map(|i| i.1)).step_by(2);
    let support = cartesian!(h_range, q_range.clone()).map(|(i, j)| 
        isize2(i, j)
    );

    Grid2::generate(support, move |idx| { 
        let i = idx.0;
        let Some(e) = info.get(&idx) else { 
            return Summand::zero()
        };
        
        let (rank, tors, indices) = e;
        let gens = grid[i].raw_gens().clone(); 
        let trans = grid[i].trans().sub(indices);
        Summand::new(gens, *rank, tors.clone(), trans)
    })
}