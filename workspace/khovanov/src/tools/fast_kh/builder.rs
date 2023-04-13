use std::collections::{VecDeque, HashSet};

use log::info;
use yui_link::{Link, Crossing, Edge};

use crate::KhComplex;

use super::complex::TngComplex;

pub struct TngComplexBuilder {
    crossings: VecDeque<Crossing>,
    complex: TngComplex
}

impl TngComplexBuilder {
    pub fn build(l: &Link) -> TngComplex { 
        info!("construct TngComplex.");

        let mut c = Self::new(l);
        c.process();

        c.complex
    }

    fn new(l: &Link) -> Self { 
        let crossings = Self::sort_crossings(l);
        let deg_shift = KhComplex::<i64>::deg_shift_for(l, false);
        let complex = TngComplex::new(deg_shift);

        Self { crossings, complex }
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
        let c = &mut self.complex;
        c.append(&x);
        
        while let Some((k, r)) = c.find_loop() { 
            let (k0, k1) = c.deloop(&k, r);
            for k in [k0, k1] { 
                if let Some((i, j)) = c.find_inv_edge(&k) { 
                    c.eliminate(&i.clone(), &j.clone());
                }
            }
        }
    }
}

#[cfg(test)]
mod tests { 
    use yui_homology::{RModStr, HomologyComputable};
    use yui_homology::test::ChainComplexValidation;
    use super::*;

    #[test]
    fn test_unknot_rm1() {
        let l = Link::from(&[[0,0,1,1]]);
        let c = TngComplexBuilder::build(&l);
        let c = c.as_generic(0, 0);

        assert_eq!(c[0].rank(), 2);
        assert_eq!(c[1].rank(), 0);

        let l = Link::from(&[[0,1,1,0]]);
        let c = TngComplexBuilder::build(&l);
        let c = c.as_generic(0, 0);

        c.check_d_all();

        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[0].rank(), 2);
    }

    #[test]
    fn test_unknot_rm2() {
        let l = Link::from(&[[1,4,2,1],[2,4,3,3]]);
        let c = TngComplexBuilder::build(&l);
        let c = c.as_generic(0, 0);

        c.check_d_all();

        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[0].rank(), 2);
        assert_eq!(c[1].rank(), 0);
    }

    #[test]
    fn test_unlink_2() {
        let pd_code = [[1,2,3,4], [3,2,1,4]];
        let l = Link::from(&pd_code);
        let c = TngComplexBuilder::build(&l);
        let c = c.as_generic(0, 0);

        c.check_d_all();

        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[0].rank(), 4);
        assert_eq!(c[1].rank(), 0);
    }

    #[test]
    fn test_hopf_link() {
        let l = Link::hopf_link();
        let c = TngComplexBuilder::build(&l);
        let c = c.as_generic(0, 0);

        c.check_d_all();

        assert_eq!(c[-2].rank(), 2);
        assert_eq!(c[-1].rank(), 0);
        assert_eq!(c[0].rank(), 2);
    }

    #[test]
    fn test_8_19() {
        let l = Link::from(&[[4,2,5,1],[8,4,9,3],[9,15,10,14],[5,13,6,12],[13,7,14,6],[11,1,12,16],[15,11,16,10],[2,8,3,7]]);
        let c = TngComplexBuilder::build(&l);
        let c = c.as_generic(0, 0);

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
}