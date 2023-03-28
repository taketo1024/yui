use std::hash::Hash;
use std::collections::HashSet;
use std::fmt::Display;
use derive_more::Display;
use itertools::Itertools;
use yui_core::Elem;
use yui_lin_comb::{FreeGen, OrdForDisplay};
use yui_utils::set;

use super::tng::{Tng, TngUpdate};

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Display)]
pub enum Dot { 
    X, Y
}

#[derive(PartialEq, Eq, Clone, Copy)]
enum End { 
    Src, Tgt
}

impl End { 
    fn is_src(&self) -> bool { 
        self == &End::Src
    }

    fn is_tgt(&self) -> bool { 
        self == &End::Tgt
    }
}

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct CobComp { 
    src: HashSet<usize>, // indices of comps in the source tangle. 
    tgt: HashSet<usize>,
    dots: (usize, usize) // nums of X and Y dots resp. 
}

impl CobComp { 
    pub fn new(src: HashSet<usize>, tgt: HashSet<usize>, dots: (usize, usize)) -> Self { 
        Self { src, tgt, dots }
    }

    pub fn cyl(i: usize, j: usize) -> Self { 
        Self::new(set![i], set![j], (0, 0))
    }

    fn contains(&self, i: usize, e: End) -> bool { 
        self.end(e).contains(&i)
    }
    
    fn connect(&mut self, other: Self) { 
        let CobComp{ src, tgt, dots } = other;

        assert_eq!(self.src.intersection(&src).count(), 1);
        assert_eq!(self.tgt.intersection(&tgt).count(), 1);

        self.src.extend(src);
        self.tgt.extend(tgt);

        self.dots.0 += dots.0;
        self.dots.1 += dots.1;
    }

    fn end(&self, e: End) -> &HashSet<usize> { 
        if e.is_src() { 
            &self.src 
        } else { 
            &self.tgt
        }
    }

    fn end_mut(&mut self, e: End) -> &mut HashSet<usize> { 
        if e.is_src() { 
            &mut self.src 
        } else { 
            &mut self.tgt
        }
    }
}

impl Display for CobComp {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let dots = vec!["X"; self.dots.0].join("") + &vec!["Y"; self.dots.1].join("");
        write!(f, "{{{:?} -> {:?}}}{}", self.src, self.tgt, dots)
    }
}

impl Hash for CobComp {
    fn hash<H: std::hash::Hasher>(&self, _state: &mut H) {
        // TODO
    }
}

#[derive(Clone, PartialEq, Eq, Hash, Debug, Default)]
pub struct Cob { 
    comps: Vec<CobComp>
}

impl Cob {
    pub fn new(comps: Vec<CobComp>) -> Self { 
        Self { comps }
    }
    
    pub fn id_for(v: &Tng) -> Self { 
        let comps = (0..v.ncomps()).map(|i| 
            CobComp::cyl(i, i)
        ).collect();
        Self::new(comps)
    }

    pub fn ncomps(&self) -> usize { 
        self.comps.len()
    }

    pub fn append_cyl(&mut self, u0: &Vec<TngUpdate>, u1: &Vec<TngUpdate>) { 
        assert_eq!(u0.len(), 2);
        assert_eq!(u1.len(), 2);

        for i in 0..2 { 
            let (u0, u1) = (&u0[i], &u1[i]);
            self.append_arc_cyl(u0, u1);
        }
    }
    
    fn append_arc_cyl(&mut self, u0: &TngUpdate, u1: &TngUpdate) { 
        assert_eq!(u0.added(), u1.added());

        self.modify_end(u0, End::Src);
        self.modify_end(u1, End::Tgt);

        let i0 = u0.index();
        let i1 = u1.index();

        if u0.added() { 
            let c = CobComp::cyl(i0, i1);
            self.comps.push(c);
        } else { 
            self.connect_if_disj(i0, End::Src);
            assert!(!self.connect_if_disj(i1, End::Tgt)) // comps should be already connected
        }
    }

    pub fn append_sdl(&mut self, u0: &Vec<TngUpdate>, u1: &Vec<TngUpdate>) { 
        assert_eq!(u0.len(), 2);
        assert_eq!(u1.len(), 2);

        // TODO
    }

    fn modify_end(&mut self, u: &TngUpdate, e: End) {
        let Some(j) = u.removed() else { return };
        let i = u.index();

        for c in self.comps.iter_mut() { 
            let end = c.end_mut(e);
            if end.contains(&i) && end.contains(&j) { 
                end.remove(&j);
            }
            *end = end.iter().map(|&k| { 
                u.apply(k)
            }).collect();
        }
    }

    fn connect_if_disj(&mut self, i: usize, e: End) -> bool { 
        let ks = (0 .. self.ncomps()).filter(|&k| 
            self.comps[k].contains(i, e)
        ).collect_vec();

        assert!(!ks.is_empty());
        assert!(ks.len() <= 2);

        if ks.len() == 2 { 
            let (k0, k1) = (ks[0], ks[1]);

            let c1 = self.comps.remove(k1);
            let c0 = &mut self.comps[k0];
            c0.connect(c1);
            
            true
        } else { 
            false
        }
    }
}

impl Display for Cob {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let cobs = self.comps.iter().map(|c| c.to_string()).join(", ");
        write!(f, "cob[{}]", cobs)
    }
}

impl OrdForDisplay for Cob {
    fn cmp_for_display(&self, _other: &Self) -> std::cmp::Ordering {
        std::cmp::Ordering::Equal
    }
}

impl Elem for Cob {
    fn set_symbol() -> String {
        "Cob".to_string()
    }
}

impl FreeGen for Cob {}

#[cfg(test)]
mod tests {
    use super::CobComp;
    use super::*;
 
    #[test]
    fn contains() { 
        let c = CobComp::new(set![0,1,2], set![3,4], (0, 0));
        assert!( c.contains(0, End::Src));
        assert!(!c.contains(3, End::Src));
        assert!( c.contains(3, End::Tgt));
        assert!(!c.contains(0, End::Tgt));
    }

    #[test]
    fn connect() { 
        let mut c0 = CobComp::new(set![0,1,2], set![10,11],    (1, 2));
        let     c1 = CobComp::new(set![2,3],   set![11,12,13], (3, 1));

        c0.connect(c1);

        assert_eq!(c0.src, set![0,1,2,3]);
        assert_eq!(c0.tgt, set![10,11,12,13]);
        assert_eq!(c0.dots, (4, 3));
    }

    #[test]
    fn modify_end() { 
        let c0 = Cob::new(vec![
            CobComp::new(set![0,1,2], set![10,11], (0, 0)),
            CobComp::new(set![3,4,5], set![12,13], (0, 0))
        ]);
        
        let mut c = c0.clone();
        let u = TngUpdate::new(1, false, None);
        c.modify_end(&u, End::Src); // nothing happens
        
        assert_eq!(c, c0);

        let mut c = c0.clone();
        let u = TngUpdate::new(1, false, Some(2));
        c.modify_end(&u, End::Src);
        
        assert_eq!(c, Cob::new(vec![
            CobComp::new(set![0,1], set![10,11], (0, 0)),
            CobComp::new(set![2,3,4], set![12,13], (0, 0))
        ]));

        let mut c = c0.clone();
        let u = TngUpdate::new(2, false, Some(3));
        c.modify_end(&u, End::Src);
        
        assert_eq!(c, Cob::new(vec![
            CobComp::new(set![0,1,2], set![10,11], (0, 0)),
            CobComp::new(set![3,2,4], set![12,13], (0, 0))
        ]));

        let mut c = c0.clone();
        let u = TngUpdate::new(10, false, Some(12));
        c.modify_end(&u, End::Tgt);
        
        assert_eq!(c, Cob::new(vec![
            CobComp::new(set![0,1,2], set![10,11], (0, 0)),
            CobComp::new(set![3,4,5], set![10,12], (0, 0))
        ]));
    }

    #[test]
    fn connect_if_disj() { 
        let c0 = Cob::new(vec![
            CobComp::new(set![0,1,2], set![10,11], (0, 0)),
            CobComp::new(set![2,3,4], set![11,12], (0, 0))
        ]);

        let mut c = c0.clone();
        let r = c.connect_if_disj(1, End::Src); // nothing happens
        
        assert_eq!(c, c0);
        assert_eq!(r, false);

        let mut c = c0.clone();
        let r = c.connect_if_disj(2, End::Src);

        assert_eq!(c, Cob::new(vec![
            CobComp::new(set![0,1,2,3,4], set![10,11,12], (0, 0)),
        ]));
        assert_eq!(r, true);

        let mut c = c0.clone();
        let r = c.connect_if_disj(11, End::Tgt);

        assert_eq!(c, Cob::new(vec![
            CobComp::new(set![0,1,2,3,4], set![10,11,12], (0, 0)),
        ]));
        assert_eq!(r, true);
    }

    #[test]
    fn append_arc_cyl() { 
        let c0 = Cob::new(vec![
            CobComp::new(set![0,1,2], set![10,11], (0, 0)),
            CobComp::new(set![2,3,4], set![11,12], (0, 0))
        ]);
        
        // when arcs are newly added.
        
        let mut c = c0.clone();
        let u0 = TngUpdate::new(5,  true, None);
        let u1 = TngUpdate::new(13, true, None);
        c.append_arc_cyl(&u0, &u1);

        assert_eq!(c, Cob::new(vec![
            CobComp::new(set![0,1,2], set![10,11], (0, 0)),
            CobComp::new(set![2,3,4], set![11,12], (0, 0)),
            CobComp::new(set![5],     set![13],    (0, 0)),            
        ]));

        // when arcs connect ends.
        
        let mut c = c0.clone();
        let u0 = TngUpdate::new(2,  false, None);
        let u1 = TngUpdate::new(11, false, None);
        c.append_arc_cyl(&u0, &u1);

        assert_eq!(c, Cob::new(vec![
            CobComp::new(set![0,1,2,3,4], set![10,11,12], (0, 0)),
        ]));
    }
}