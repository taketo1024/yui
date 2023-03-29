use std::hash::Hash;
use std::collections::HashSet;
use std::fmt::Display;
use derive_more::Display;
use itertools::Itertools;
use yui_core::Elem;
use yui_lin_comb::{FreeGen, OrdForDisplay};
use yui_polynomial::Mono2;
use yui_utils::set;

use super::tng::{Tng, TngUpdate};

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Display)]
pub enum Dot { 
    None, X, Y
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

    pub fn zero() -> Self { // sphere without dots
        Self::new(set![], set![], (0, 0))
    }

    pub fn cyl(i0: usize, i1: usize) -> Self { 
        Self::new(
            set![i0], 
            set![i1], 
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

    pub fn contains(&self, i: usize, e: End) -> bool { 
        self.end(e).contains(&i)
    }

    pub fn is_connectable(&self, other: &Self) -> bool { 
        self.src.intersection(&other.src).next().is_some() &&
        self.tgt.intersection(&other.tgt).next().is_some()
    }
    
    pub fn connect(&mut self, other: Self) { 
        let CobComp{ src, tgt, dots } = other;

        self.src.extend(src);
        self.tgt.extend(tgt);

        self.dots.0 += dots.0;
        self.dots.1 += dots.1;
    }

    pub fn has_dots(&self) -> bool { 
        self.dots.0 > 0 || self.dots.1 > 0
    }

    pub fn is_closed(&self) -> bool { // sphere possibly with dots.
        self.src.is_empty() && self.tgt.is_empty()
    }

    pub fn is_zero(&self) -> bool { 
        self.is_closed() && self.dots == (0, 0)
    }

    pub fn is_removable(&self) -> bool { 
        self.is_closed() && (self.dots == (1, 0) || self.dots == (0, 1))
    }

    pub fn cup(&mut self, i0: usize) {
        assert!( self.src.remove(&i0) )
    }

    pub fn cap(&mut self, i1: usize) {
        assert!( self.tgt.remove(&i1) )
    }

    pub fn dot(&mut self, dot: Dot) { 
        match dot { 
            Dot::X => self.dots.0 += 1,
            Dot::Y => self.dots.1 += 1,
            _      => ()
        }
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
        fn set(s: &HashSet<usize>) -> String { 
            if s.is_empty() { 
                "∅".to_string()
            } else { 
                format!("{{{}}}", s.iter().sorted().join(","))
            }
        }

        let s = match (self.src.len(), self.tgt.len()) { 
            (0, 0) => "S",
            (0, 1) => "ι",
            (1, 0) => "ε",
            (1, 1) => "id",
            (2, 1) => "m",
            (1, 2) => "Δ",
            _      => "cob"
        };

        let dots = if self.has_dots() { 
            let (p, q) = self.dots;
            format!("; {}", Mono2::<'X','Y', _>::new(p, q).to_string())
        } else { 
            String::new()
        };
        
        write!(f, "{s}({} -> {}{})", set(&self.src), set(&self.tgt), dots)
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

    pub fn is_zero(&self) -> bool { 
        self.comps.iter().find(|c| c.is_zero()).is_some()
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

    pub fn insert(&mut self, mut c: CobComp) { // horizontal composition
        let mut i = 0;

        while i < self.comps.len() { 
            if c.is_connectable(&self.comps[i]) { 
                let c2 = self.comps.remove(i);
                c.connect(c2);
            } else { 
                i += 1;
            }
        }
        
        self.comps.push(c);

        #[cfg(debug_assertions)] {
            self.comps.sort_by_key(|c| c.src.iter().min().cloned().unwrap_or(0));
        }
    }

    pub fn cup(&mut self, r: usize, x: Dot) {
        self.cup_or_cap(r, x, End::Src)
    }

    pub fn cap(&mut self, r: usize, x: Dot) {
        self.cup_or_cap(r, x, End::Tgt)
    }

    pub fn cup_or_cap(&mut self, r: usize, x: Dot, e: End) {
        let (i, c) = self.comp_containing(r, e);

        if e.is_src() { 
            c.cup(r);
        } else {
            c.cap(r);
        }
        c.dot(x);

        if c.is_zero() { 
            self.comps.clear();
            self.comps.push(CobComp::zero())
        } else if c.is_removable() { 
            self.comps.remove(i);
        }

        // reindex
        for c in self.comps.iter_mut() { 
            let end = c.end_mut(e);
            *end = end.iter().map(|&k| { 
                if k > r { 
                    k - 1
                } else { 
                    k
                }
            }).collect();
        }
    }

    fn comp_containing(&mut self, r: usize, e: End) -> (usize, &mut CobComp) { 
        self.comps.iter_mut().enumerate().find(|(_, c)| 
            c.contains(r, e)
        ).unwrap()
    }
}

impl Display for Cob {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.comps.is_empty() { 
            write!(f, "1")
        } else { 
            let cobs = self.comps.iter().map(|c| c.to_string()).join(" ⊔ ");
            write!(f, "[{}]", cobs)
        }
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
    fn is_connectable() { 
        let c0 = CobComp::new(set![0,1,2], set![10,11],    (0, 0));
        let c1 = CobComp::new(set![2,3],   set![11,12,13], (0, 0));
        let c2 = CobComp::new(set![3,4,5], set![13,14], (0, 0));

        assert!( c0.is_connectable(&c1));
        assert!( c1.is_connectable(&c2));
        assert!(!c0.is_connectable(&c2));
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
    fn apply_update() { 
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
    fn insert() { 
        let mut c = Cob::new(vec![]);
        c.insert(CobComp::new(set![0,1,2], set![10,11], (0, 0)));
        c.insert(CobComp::new(set![3,4,5], set![12,13], (0, 0)));

        assert_eq!(c.comps.len(), 2);
        
        c.insert(CobComp::new(set![2,3], set![11,12], (0, 0)));

        assert_eq!(c.comps.len(), 1);
        assert_eq!(c, Cob::new(vec![
            CobComp::new(set![0,1,2,3,4,5], set![10,11,12,13], (0, 0)),
        ]));
    }
}