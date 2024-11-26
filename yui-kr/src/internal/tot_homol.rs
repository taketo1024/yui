use std::ops::{Index, RangeInclusive};
use std::sync::Arc;
use cartesian::cartesian;
use delegate::delegate;

use itertools::Itertools;
use log::info;
use yui::{EucRing, EucRingOps};
use yui_homology::{isize2, ChainComplexTrait, ComputeHomology, DisplayTable, GenericHomology2, GenericSummand, Grid2, GridIter, GridTrait};

use super::tot_cpx::KRTotComplex;
use super::base::extend_ends_bounded;
use super::data::KRCubeData;

#[derive(Default)]
pub struct KRTotHomol<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    q: isize,
    inner: GenericHomology2<R>
} 

impl<R> KRTotHomol<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    pub fn new(data: Arc<KRCubeData<R>>, q: isize) -> Self { 
        let n = data.dim() as isize;
        Self::new_restr(data, q, (0..=n, 0..=n))
    }

    pub fn new_restr(data: Arc<KRCubeData<R>>, q: isize, range: (RangeInclusive<isize>, RangeInclusive<isize>)) -> Self { 
        info!("H (q: {}, h: {:?}, v: {:?})..", q, range.0, range.1);

        let support = cartesian!(
            range.0.clone(), 
            range.1.clone()
        ).map(isize2::from);

        let complex = Self::make_cpx(data, q, range.clone());
        let complex = complex.reduced();

        let inner = Grid2::generate(
            support,
            |idx| complex.compute_homology_at(idx, false)
        );

        info!("H (q: {}, h: {:?}, v: {:?})\n{}", q, range.0, range.1, inner.display_table("h", "v"));

        Self { q, inner }
    }

    pub fn try_partial(data: Arc<KRCubeData<R>>, q: isize, range: (RangeInclusive<isize>, RangeInclusive<isize>), size_limit: usize) -> Grid2<Option<GenericSummand<isize2, R>>> {
        info!("H (q: {}, h: {:?}, v: {:?})..", q, range.0, range.1);

        let support = cartesian!(
            range.0.clone(), 
            range.1.clone()
        ).map(isize2::from);

        let complex = Self::make_cpx(data, q, range.clone());
        let complex = complex.reduced_with_limit(size_limit);

        let d_deg = complex.d_deg();
        let grid = Grid2::generate(
            support,
            |idx| (complex.rank(idx - d_deg) <= size_limit 
                && complex.rank(idx)         <= size_limit 
                && complex.rank(idx + d_deg) <= size_limit).then(|| 
                    complex.compute_homology_at(idx, false)
                )
        );

        info!("H (q: {}, h: {:?}, v: {:?})\n{}", q, range.0, range.1, 
            display_table(&grid, "h", "v", |h| h.as_ref().map(|h| h.to_string()).unwrap_or("?".to_string()))
        );

        grid
    }

    fn make_cpx(data: Arc<KRCubeData<R>>, q: isize, range: (RangeInclusive<isize>, RangeInclusive<isize>)) -> KRTotComplex<R> { 
        let n = data.dim() as isize;
        let c_range = (
            range.0.clone(), 
            extend_ends_bounded(range.1.clone(), 1, 0..=n)
        );
        
        KRTotComplex::new_restr(
            data.clone(), 
            q, 
            c_range
        )
    }

    pub fn q_deg(&self) -> isize { 
        self.q
    }
}

impl<R> GridTrait<isize2> for KRTotHomol<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    type Support = GridIter<isize2>;
    type Item = GenericSummand<isize2, R>;

    delegate! { 
        to self.inner {
            fn support(&self) -> Self::Support;
            fn is_supported(&self, i: isize2) -> bool;
            fn get(&self, idx: isize2) -> &Self::Item;
            fn get_default(&self) -> &Self::Item;
        }
    }
}

impl<R> Index<(isize, isize)> for KRTotHomol<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    type Output = GenericSummand<isize2, R>;
    fn index(&self, index: (isize, isize)) -> &Self::Output { 
        self.get(index.into())
    }
}

// TODO to be removed
fn display_table<E, F>(grid: &Grid2<E>, label0: &str, label1: &str, f: F) -> String
where F: Fn(&E) -> String {
    use yui::util::format::table;

    let head = format!("{}\\{}", label1, label0);
    let cols = grid.support().map(|e| e.0).unique().sorted();
    let rows = grid.support().map(|e| e.1).unique().sorted().rev();

    let str = table(head, rows, cols, |&j, &i| {
        f(grid.get(isize2(i, j)))
    });

    str
}

#[cfg(test)]
mod tests { 
    use yui_homology::SummandTrait;
    use yui_link::Link;
    use yui::Ratio;
    use super::*;

    type R = Ratio<i64>;

    #[test]
    fn rank() { 
        let l = Link::from_pd_code([[1,4,2,5],[5,2,6,3],[3,6,4,1]]); // trefoil
        let q = -4;
        let data = Arc::new( KRCubeData::<R>::new(&l) );
        let h = KRTotHomol::new(data, q);

        assert_eq!(h[(0, 0)].rank(), 0);
        assert_eq!(h[(0, 1)].rank(), 0);
        assert_eq!(h[(0, 2)].rank(), 0);
        assert_eq!(h[(0, 3)].rank(), 0);
        assert_eq!(h[(1, 0)].rank(), 0);
        assert_eq!(h[(1, 1)].rank(), 0);
        assert_eq!(h[(1, 2)].rank(), 0);
        assert_eq!(h[(1, 3)].rank(), 0);
        assert_eq!(h[(2, 0)].rank(), 0);
        assert_eq!(h[(2, 1)].rank(), 0);
        assert_eq!(h[(2, 2)].rank(), 0);
        assert_eq!(h[(2, 3)].rank(), 1);
        assert_eq!(h[(3, 0)].rank(), 0);
        assert_eq!(h[(3, 1)].rank(), 1);
        assert_eq!(h[(3, 2)].rank(), 0);
        assert_eq!(h[(3, 3)].rank(), 0);
     }
}