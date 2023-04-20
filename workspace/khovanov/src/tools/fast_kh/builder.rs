use std::collections::{VecDeque, HashSet};

use log::info;
use yui_link::{Link, Crossing, Edge};

use crate::{KhComplex, KhEnhState, KhLabel};
use crate::tools::fast_kh::cob::{CobComp, Dot};
use crate::tools::fast_kh::tng::{TngComp};

use super::cob::Cob;
use super::complex::TngComplex;
use super::elem::TngElem;

pub struct TngComplexBuilder {
    crossings: VecDeque<Crossing>,
    complex: TngComplex,
    canon_cycles: Vec<TngElem>
}

impl TngComplexBuilder {
    pub fn build(l: &Link, with_canon_cycle: bool) -> Self { 
        info!("construct TngComplex.");

        let mut c = Self::new(l);
        
        if with_canon_cycle { 
            assert_eq!(l.components().len(), 1);
            c.make_canon_cycles();
        }
        
        c.process();

        c
    }

    fn new(l: &Link) -> Self { 
        let crossings = Self::sort_crossings(l);
        let deg_shift = KhComplex::<i64>::deg_shift_for(l, false);
        let complex = TngComplex::new(deg_shift);

        Self { crossings, complex, canon_cycles: vec![] }
    }

    fn sort_crossings(l: &Link) -> VecDeque<Crossing> { 
        let mut remain = l.data().clone();
        let mut endpts = HashSet::new();
        let mut res = VecDeque::new();

        fn take_best(remain: &mut Vec<Crossing>, endpts: &mut HashSet<Edge>) -> Option<Crossing> { 
            if remain.is_empty() { 
                return None 
            }

            let mut cand_i = 0;
            let mut cand_c = 0;

            for (i, x) in remain.iter().enumerate() { 
                let c = x.edges().iter().filter(|e| 
                    endpts.contains(e)
                ).count();

                if c == 4 { 
                    let x = remain.remove(i);
                    return Some(x);
                } else if c > cand_c { 
                    cand_i = i;
                    cand_c = c;
                }
            }

            let x = remain.remove(cand_i);
            for e in x.edges() { 
                endpts.insert(*e);
            }
            
            Some(x)
        }

        while let Some(x) = take_best(&mut remain, &mut endpts) { 
            res.push_back(x);
        }

        res
    }

    fn process(&mut self) {
        while let Some(x) = self.crossings.pop_front() { 
            self.proceed_each(x);
        }
        // info!("result: #v = {}", c.vertices.len());
    }

    fn proceed_each(&mut self, x: Crossing) { 
        self.complex.append(&x);

        while let Some((k, r)) = self.complex.find_loop() { 
            self.deloop(&k, r);
        }
    }

    fn deloop(&mut self, k: &KhEnhState, r: usize) {
        let (c, k0, k1) = self.complex.deloop(k, r);

        for e in self.canon_cycles.iter_mut() { 
            e.deloop(k, &c);
        }

        for k in [k0, k1] { 
            self.eliminate(&k)
        }
    }

    fn eliminate(&mut self, k: &KhEnhState) {
        if let Some((i, j)) = self.complex.find_inv_edge(k) { 
            self.complex.eliminate(&i, &j);
        }
    }

    fn make_canon_cycles(&mut self) { 
        use crate::derived::canon_cycle::ColoredSeifertCircles;

        let l = Link::new(self.crossings.iter().cloned().collect());
        assert_eq!(l.components().len(), 1);

        let s = l.ori_pres_state();
        let circles = l.colored_seifert_circles(true);

        let t = KhEnhState::new(s, KhLabel::empty());
        let f = Cob::new(
            circles.iter().map(|(circ, col)| { 
                let t = TngComp::from(circ);
                let mut cup = CobComp::cup(t);
                let dot = if col.is_a() { Dot::X } else { Dot::Y };
                cup.add_dot(dot);
                cup
            }).collect()
        );
        let e = TngElem::init(t, f);

        info!("canon-cycle: {}", e);
        
        self.canon_cycles = vec![e];
    }
}

#[cfg(test)]
mod tests { 
    use yui_homology::{RModStr, HomologyComputable};
    use yui_homology::test::ChainComplexValidation;
    use super::*;

    #[test]
    fn test_unknot_rm1() {
        let l = Link::from_pd_code([[0,0,1,1]]);
        let b = TngComplexBuilder::build(&l, false);
        let c = b.complex.as_generic(0, 0);

        assert_eq!(c[0].rank(), 2);
        assert_eq!(c[1].rank(), 0);

        let l = Link::from_pd_code([[0,1,1,0]]);
        let b = TngComplexBuilder::build(&l, false);
        let c = b.complex.as_generic(0, 0);

        c.check_d_all();

        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[0].rank(), 2);
    }

    #[test]
    fn test_unknot_rm2() {
        let l = Link::from_pd_code([[1,4,2,1],[2,4,3,3]]);
        let b = TngComplexBuilder::build(&l, false);
        let c = b.complex.as_generic(0, 0);

        c.check_d_all();

        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[0].rank(), 2);
        assert_eq!(c[1].rank(), 0);
    }

    #[test]
    fn test_unlink_2() {
        let pd_code = [[1,2,3,4], [3,2,1,4]];
        let l = Link::from_pd_code(pd_code);
        let b = TngComplexBuilder::build(&l, false);
        let c = b.complex.as_generic(0, 0);

        c.check_d_all();

        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[0].rank(), 4);
        assert_eq!(c[1].rank(), 0);
    }

    #[test]
    fn test_hopf_link() {
        let l = Link::hopf_link();
        let b = TngComplexBuilder::build(&l, false);
        let c = b.complex.as_generic(0, 0);

        c.check_d_all();

        assert_eq!(c[-2].rank(), 2);
        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[0].rank(), 2);
    }

    #[test]
    fn test_8_19() {
        let l = Link::from_pd_code([[4,2,5,1],[8,4,9,3],[9,15,10,14],[5,13,6,12],[13,7,14,6],[11,1,12,16],[15,11,16,10],[2,8,3,7]]);
        let b = TngComplexBuilder::build(&l, false);
        let c = b.complex.as_generic(0, 0);

        c.check_d_all();

        let h = c.homology();

        for i in [1,6,7,8] {
            assert_eq!(h[i].rank(), 0);
            assert_eq!(h[i].is_free(), true);
        }

        for i in [0,4,5] {
            assert_eq!(h[i].rank(), 2);
            assert_eq!(h[i].is_free(), true);
        }

        assert_eq!(h[2].rank(), 1);
        assert_eq!(h[2].is_free(), true);

        assert_eq!(h[3].rank(), 1);
        assert_eq!(h[3].tors(), &vec![2]);
    }

    #[test]
    fn canon_cycle_trefoil() { 
        let l = Link::trefoil();
        let b = TngComplexBuilder::build(&l, true);
        // println!("{}", b.canon_cycles[0]);
    }
}