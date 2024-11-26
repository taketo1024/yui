use std::ops::{Index, RangeInclusive};
use std::sync::Arc;

use delegate::delegate;
use yui::{EucRing, EucRingOps};
use yui_homology::{isize2, isize3, GenericSummand, Grid1, GridIter, GridTrait, SummandTrait};
use yui_link::{Link, Braid};

use crate::KRHomologyStr;
use crate::internal::data::KRCubeData;
use crate::internal::tot_homol::KRTotHomol;

pub struct KRHomology<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    data: Arc<KRCubeData<R>>,
    q_slices: Grid1<KRTotHomol<R>>
}

impl<R> KRHomology<R> 
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    pub fn from_braid(b: &Braid) -> Self { 
        let l = b.closure();
        Self::from_link(&l)
    }

    pub fn from_link(link: &Link) -> Self { 
        assert_eq!(link.components().len(), 1, "only knots are supported.");
        let data = Arc::new(KRCubeData::new(link));
        Self::new(data)
    }

    pub fn new(data: Arc<KRCubeData<R>>) -> Self { 
        let q_range = data.q_range();
        let q_slices = Grid1::generate(q_range, |q| {
            let range = data.hv_range(q);
            KRTotHomol::new_restr(data.clone(), q, range)
        });
        Self { data, q_slices }
    }

    delegate! { 
        to self.data { 
            pub fn i_range(&self) -> RangeInclusive<isize>;
            pub fn j_range(&self) -> RangeInclusive<isize>;
            pub fn k_range(&self) -> RangeInclusive<isize>;
            pub fn q_range(&self) -> RangeInclusive<isize>;
        }
    }

    pub fn structure(&self) -> KRHomologyStr { 
        self.support().filter_map(|idx| {
            let h = self.get(idx);
            let r = h.rank();

            if r > 0 { 
                Some(((idx.0, idx.1, idx.2), r))
            } else { 
                None
            }
        }).collect()
    }

    pub fn display_table(&self) -> String { 
        self.structure().qpoly_table()
    }

    pub fn print_table(&self) { 
        println!("{}", self.display_table())
    }
}

impl<R> GridTrait<isize3> for KRHomology<R> 
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    type Support = GridIter<isize3>;
    type Item = GenericSummand<isize2, R>;

    delegate! { 
        to self.data { 
            fn support(&self) -> Self::Support;
            fn is_supported(&self, idx: isize3) -> bool;
        }
    }

    fn get(&self, idx: isize3) -> &Self::Item {
        if let Some(isize3(q, h, v)) = self.data.to_inner_grad(idx) { 
            &self.q_slices[q][(h, v)]
        } else { 
            self.get_default()
        }
    }
    
    fn get_default(&self) -> &Self::Item {
        self.q_slices.get_default().get(isize2(0, 0))
    }
}

impl<R> Index<(isize, isize, isize)> for KRHomology<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    type Output = GenericSummand<isize2, R>;
    fn index(&self, index: (isize, isize, isize)) -> &Self::Output {
        self.get(index.into())
    }
}

#[cfg(test)]
mod tests { 
    use yui_link::Braid;
    use yui::Ratio;
    use yui::util::macros::hashmap;
    use super::*;

    type R = Ratio<i64>;

    #[test]
    fn trivial() { 
        let b = Braid::from([1]);
        let h = KRHomology::<R>::from_braid(&b);

        assert_eq!(h[(0,0,0)].rank(), 1);
        assert_eq!(h[(1,0,0)].rank(), 0);

        let s = h.structure();

        assert_eq!(s.total_rank(), 1);
        assert_eq!(s, KRHomologyStr::from(hashmap!{ 
            (0,0,0) => 1
        }));
    }

    #[test]
    fn b3_1() { 
        let b = Braid::from([1,1,1]);
        let h = KRHomology::<R>::from_braid(&b);
        let s = h.structure();

        assert_eq!(s.total_rank(), 3);
        assert_eq!(s, KRHomologyStr::from(hashmap!{ 
            (0,4,-2) => 1,
            (-2,2,2) => 1,
            (2,2,-2) => 1
        }));
    }
    
    #[test]
    fn b4_1() { 
        let b = Braid::from([1,-2,1,-2]);
        let h = KRHomology::<R>::from_braid(&b);
        let s = h.structure();

        assert_eq!(s.total_rank(), 5);
        assert_eq!(s, KRHomologyStr::from(hashmap!{ 
            (0,-2,2) => 1,
            (-2,0,2) => 1,
            (0,2,-2) => 1,
            (0,0,0)  => 1,
            (2,0,-2) => 1
        }));
    }

    #[test]
    fn b5_1() { 
        let b = Braid::from([1,1,1,1,1]);
        let h = KRHomology::<R>::from_braid(&b);
        let s = h.structure();

        assert_eq!(s.total_rank(), 5);
        assert_eq!(s, KRHomologyStr::from(hashmap!{             
            (0,4,0)  => 1,
            (-2,6,0) => 1,
            (-4,4,4) => 1,
            (4,4,-4) => 1,
            (2,6,-4) => 1
        }));
    }

    #[test]
    fn b5_2() { 
        let b = Braid::from([1,1,1,2,-1,2]);
        let h = KRHomology::<R>::from_braid(&b);
        let s = h.structure();

        assert_eq!(s.total_rank(), 7);
        assert_eq!(s, KRHomologyStr::from(hashmap!{ 
            (2,4,-4) => 1,
            (2,2,-2) => 1,
            (0,4,-2) => 1,
            (-2,4,0) => 1,
            (0,2,0)  => 1,
            (0,6,-4) => 1,
            (-2,2,2) => 1
        }));
    }

    #[test]
    fn b6_1() { 
        let b = Braid::from([1,1,2,-1,-3,2,-3]);
        let h = KRHomology::<R>::from_braid(&b);
        let s = h.structure();

        assert_eq!(s.total_rank(), 9);
        assert_eq!(s, KRHomologyStr::from(hashmap!{ 
            (0,-2,2) => 1,
            (0,0,0)  => 2,
            (2,0,-2) => 1,
            (2,2,-4) => 1,
            (0,4,-4) => 1,
            (0,2,-2) => 1,
            (-2,2,0) => 1,
            (-2,0,2) => 1    
        }));
    }

    #[test]
    fn b6_2() { 
        let b = Braid::from([1,1,1,-2,1,-2]);
        let h = KRHomology::<R>::from_braid(&b);
        let s = h.structure();

        assert_eq!(s.total_rank(), 11);
        assert_eq!(s, KRHomologyStr::from(hashmap!{ 
            (0,4,-2) => 1,
            (2,4,-4) => 1,
            (-2,4,0) => 1,
            (0,2,0)  => 2,
            (2,2,-2) => 1,
            (2,0,0)  => 1,
            (-2,2,2) => 1,
            (-2,0,4) => 1,
            (4,2,-4) => 1,
            (-4,2,4) => 1
        }));
    }

    #[test]
    fn b6_3() { 
        let b = Braid::from([1,1,-2,1,-2,-2]);
        let h = KRHomology::<R>::from_braid(&b);
        let s = h.structure();

        assert_eq!(s.total_rank(), 13);
        assert_eq!(s, KRHomologyStr::from(hashmap!{ 
            (4,0,-4)  => 1,
            (2,-2,0)  => 1,
            (0,0,0)   => 3,
            (2,2,-4)  => 1,
            (-2,2,0)  => 1,
            (0,2,-2)  => 1,
            (0,-2,2)  => 1,
            (-2,0,2)  => 1,
            (-4,0,4)  => 1,
            (2,0,-2)  => 1,
            (-2,-2,4) => 1
        }));
    }
}