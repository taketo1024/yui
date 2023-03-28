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
pub enum End { 
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
        Self::new(
            set![i], 
            set![j], 
            (0, 0)
        )
    }

    pub fn sdl(r0: (usize, usize), r1: (usize, usize)) -> Self { 
        Self::new(
            set![r0.0, r0.1], 
            set![r1.0, r1.1], 
            (0, 0)
        )
    }

    fn contains(&self, i: usize, e: End) -> bool { 
        self.end(e).contains(&i)
    }
    
    fn connect(&mut self, other: Self) { 
        let CobComp{ src, tgt, dots } = other;

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
        write!(f, "cob({:?} -> {:?}){}", self.src, self.tgt, dots)
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

    pub fn append(&mut self, c: CobComp) { 
        self.comps.push(c);
    }

    pub fn apply_update(&mut self, u: &TngUpdate, e: End) {
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

    pub fn insert_arc_cyl(&mut self, i0: usize, i1: usize) { 
        self.connect_if_disj(i0, End::Src);
        
        debug_assert!(
            !self.connect_if_disj(i1, End::Tgt)
        )
    }

    pub fn insert_sdl(&mut self, r0: (usize, usize), r1: (usize, usize)) { 
        let (i0, j0) = r0;
        let (i1, j1) = r1;
        
        // MEMO: (at most) one of the arcs might not be included yet.

        let complement = [ 
            self.complement(i0, j0, End::Src),
            self.complement(j0, i0, End::Src),
            self.complement(i1, j1, End::Tgt),
            self.complement(j1, i1, End::Tgt),
        ];

        assert!(complement.into_iter().filter(|&b| b).count() <= 1);

        self.connect_if_disj(i0, End::Src);
        self.connect_if_disj(j0, End::Src);
        self.connect_if_disj(i1, End::Tgt);

        debug_assert!(
            !self.connect_if_disj(j1, End::Tgt)
        )
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

    fn complement(&mut self, i: usize, j: usize, e: End) -> bool { 
        if self.comp_containing(i, e).is_some() { 
            return false 
        }

        let Some(c) = self.comp_containing(j, e) else { 
            panic!() 
        };
        c.end_mut(e).insert(i);
        true
    }

    fn comp_containing(&mut self, i: usize, e: End) -> Option<&mut CobComp> { 
        self.comps.iter_mut().find(|c| c.contains(i, e))
    }
}

impl Display for Cob {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let cobs = self.comps.iter().map(|c| c.to_string()).join(" + ");
        write!(f, "[{}]", cobs)
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
        c.apply_update(&u, End::Src); // nothing happens
        
        assert_eq!(c, c0);

        let mut c = c0.clone();
        let u = TngUpdate::new(1, false, Some(2));
        c.apply_update(&u, End::Src);
        
        assert_eq!(c, Cob::new(vec![
            CobComp::new(set![0,1], set![10,11], (0, 0)),
            CobComp::new(set![2,3,4], set![12,13], (0, 0))
        ]));

        let mut c = c0.clone();
        let u = TngUpdate::new(2, false, Some(3));
        c.apply_update(&u, End::Src);
        
        assert_eq!(c, Cob::new(vec![
            CobComp::new(set![0,1,2], set![10,11], (0, 0)),
            CobComp::new(set![3,2,4], set![12,13], (0, 0))
        ]));

        let mut c = c0.clone();
        let u = TngUpdate::new(10, false, Some(12));
        c.apply_update(&u, End::Tgt);
        
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
    fn insert_arc_cyl() { 
        let mut c = Cob::new(vec![
            CobComp::new(set![0,1,2], set![10,11], (0, 0)),
            CobComp::new(set![2,3,4], set![11,12], (0, 0))
        ]);
        
        c.insert_arc_cyl(2, 11);

        assert_eq!(c, Cob::new(vec![
            CobComp::new(set![0,1,2,3,4], set![10,11,12], (0, 0)),
        ]));
    }
}