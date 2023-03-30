use std::hash::Hash;
use std::collections::HashSet;
use std::fmt::Display;
use derive_more::Display;
use itertools::Itertools;
use yui_core::{Elem, Ring, RingOps};
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
    genus: usize,
    dots: (usize, usize) // nums of X and Y dots resp. 
}

impl CobComp { 
    pub fn new(src: HashSet<usize>, tgt: HashSet<usize>, genus: usize, dots: (usize, usize)) -> Self { 
        Self { src, tgt, genus, dots }
    }

    pub fn new_elem(src: HashSet<usize>, tgt: HashSet<usize>) -> Self { 
        Self::new(src, tgt, 0, (0, 0))
    }

    pub fn cyl(i0: usize, i1: usize) -> Self { 
        Self::new_elem(
            set![i0], 
            set![i1],
        )
    }

    pub fn sdl(r0: (usize, usize), r1: (usize, usize)) -> Self { 
        Self::new_elem(
            set![r0.0, r0.1], 
            set![r1.0, r1.1], 
        )
    }

    pub fn cup(i0: usize) -> Self { 
        Self::new_elem(
            set![i0],
            set![]
        )
    }

    pub fn cap(i1: usize) -> Self { 
        Self::new_elem(
            set![],
            set![i1]
        )
    }

    pub fn genus(&self) -> usize { 
        self.genus
    }

    pub fn contains(&self, i: usize, e: End) -> bool { 
        self.end(e).contains(&i)
    }

    pub fn is_connectable(&self, other: &Self) -> bool { 
        self.src.intersection(&other.src).next().is_some() &&
        self.tgt.intersection(&other.tgt).next().is_some()
    }
    
    pub fn connect(&mut self, other: Self) { 
        let CobComp{ src, tgt, genus, dots } = other;

        self.src.extend(src);
        self.tgt.extend(tgt);

        // TODO must compute genus properly!

        self.dots.0 += dots.0;
        self.dots.1 += dots.1;
    }

    pub fn has_dots(&self) -> bool { 
        self.dots.0 > 0 || 
        self.dots.1 > 0
    }

    pub fn is_closed(&self) -> bool {
        self.src.is_empty() && 
        self.tgt.is_empty()
    }

    pub fn is_sph(&self) -> bool {
        self.is_closed() && 
        self.genus == 0
    }

    pub fn is_cyl(&self) -> bool { // (arc or circle) × I
        self.src.len() == 1 && 
        self.tgt.len() == 1 && 
        self.genus == 0
    }

    pub fn is_zero(&self) -> bool { 
        self.is_sph() && 
        (self.dots == (0, 0) || // ε.ι = 0,
         self.dots == (1, 1))   // ε.XY.ι = ε.T.ι = 0.
    }

    pub fn is_removable(&self) -> bool { 
        self.is_sph() && 
        (self.dots == (1, 0) || // ε.X.ι = 1,
         self.dots == (0, 1))   // ε.Y.ι = 1.
    }

    pub fn is_invertible(&self) -> bool { 
        self.is_cyl() && 
        self.dots == (0, 0)
    }

    pub fn inv(&self) -> Option<Self> { 
        if self.is_invertible() { 
            let &i0 = self.src.iter().next().unwrap();
            let &i1 = self.tgt.iter().next().unwrap();
            let inv = Self::new(set![i1], set![i0], 0, (0, 0));
            Some(inv)
        } else {
            None
        }
    }

    pub fn cap_off(&mut self, i: usize, e: End) {
        assert!( self.end_mut(e).remove(&i) )
    }

    pub fn add_dot(&mut self, dot: Dot) { 
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

    pub fn eval<R>(&self, h: &R, t: &R) -> R
    where R: Ring, for<'x> &'x R: RingOps<R> {
        assert!(self.is_closed());

        // TODO must `neck cut` until g = 0!

        fn eval<R>(x: usize, y: usize, h: &R, t: &R) -> R
        where R: Ring, for<'x> &'x R: RingOps<R> { 
            match (x, y) { 
                (0, 0) | (1, 1)  => R::zero(),
                (1, 0) | (0, 1)  => R::one(),
                (x, y) if x >= 2 => // X^2 = hX + t
                    h * eval(x-1, y, h, t) + 
                    t * eval(x-2, y, h, t),
                (x, y) if y >= 2 => // Y^2 = -hY + t
                    -h * eval(x, y-1, h, t) + 
                     t * eval(x, y-2, h, t),
                _ => panic!()
            }
        }

        let (x, y) = self.dots;
        eval(x, y, h, t)
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

    pub fn is_closed(&self) -> bool { 
        self.comps.iter().all(|c| c.is_closed())
    }

    pub fn is_invertible(&self) -> bool { 
        self.comps.iter().all(|c| c.is_invertible())
    }

    pub fn inv(&self) -> Option<Self> { 
        if self.is_invertible() { 
            let comps = self.comps.iter().map(|c| c.inv().unwrap()).collect();
            let inv = Self::new(comps);
            Some(inv)
        } else { 
            None
        }
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
                u.reindex(k)
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

    pub fn cap_off(&mut self, r: usize, x: Dot, e: End) {
        let (i, c) = self.comp_containing(r, e);

        c.cap_off(r, e);
        c.add_dot(x);

        if c.is_removable() { 
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

    pub fn eval<R>(&self, h: &R, t: &R) -> R
    where R: Ring, for<'x> &'x R: RingOps<R> {
        self.comps.iter().map(|c| 
            c.eval(h, t)
        ).product()
    }
}

impl From<CobComp> for Cob {
    fn from(c: CobComp) -> Self {
        Self::new(vec![c])
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
        let c = CobComp::new(set![0,1,2], set![3,4], 0, (0, 0));
        assert!( c.contains(0, End::Src));
        assert!(!c.contains(3, End::Src));
        assert!( c.contains(3, End::Tgt));
        assert!(!c.contains(0, End::Tgt));
    }

    #[test]
    fn is_connectable() { 
        let c0 = CobComp::new(set![0,1,2], set![10,11],    0, (0, 0));
        let c1 = CobComp::new(set![2,3],   set![11,12,13], 0, (0, 0));
        let c2 = CobComp::new(set![3,4,5], set![13,14],    0, (0, 0));

        assert!( c0.is_connectable(&c1));
        assert!( c1.is_connectable(&c2));
        assert!(!c0.is_connectable(&c2));
    }

    #[test]
    fn connect() { 
        let mut c0 = CobComp::new(set![0,1,2], set![10,11],    0, (1, 2));
        let     c1 = CobComp::new(set![2,3],   set![11,12,13], 0, (3, 1));

        c0.connect(c1);

        assert_eq!(c0.src, set![0,1,2,3]);
        assert_eq!(c0.tgt, set![10,11,12,13]);
        assert_eq!(c0.dots, (4, 3));
    }

    #[test]
    fn apply_update() { 
        let c0 = Cob::new(vec![
            CobComp::new(set![0,1,2], set![10,11], 0, (0, 0)),
            CobComp::new(set![3,4,5], set![12,13], 0, (0, 0))
        ]);
        
        let mut c = c0.clone();
        let u = TngUpdate::new(1, None);
        c.apply_update(&u, End::Src); // nothing happens
        
        assert_eq!(c, c0);

        let mut c = c0.clone();
        let u = TngUpdate::new(1, Some(2));
        c.apply_update(&u, End::Src);
        
        assert_eq!(c, Cob::new(vec![
            CobComp::new(set![0,1],   set![10,11], 0, (0, 0)),
            CobComp::new(set![2,3,4], set![12,13], 0, (0, 0))
        ]));

        let mut c = c0.clone();
        let u = TngUpdate::new(2, Some(3));
        c.apply_update(&u, End::Src);
        
        assert_eq!(c, Cob::new(vec![
            CobComp::new(set![0,1,2], set![10,11], 0, (0, 0)),
            CobComp::new(set![3,2,4], set![12,13], 0, (0, 0))
        ]));

        let mut c = c0.clone();
        let u = TngUpdate::new(10, Some(12));
        c.apply_update(&u, End::Tgt);
        
        assert_eq!(c, Cob::new(vec![
            CobComp::new(set![0,1,2], set![10,11], 0, (0, 0)),
            CobComp::new(set![3,4,5], set![10,12], 0, (0, 0))
        ]));
    }

    #[test]
    fn insert() { 
        let mut c = Cob::new(vec![]);
        c.insert(CobComp::new(set![0,1,2], set![10,11], 0, (0, 0)));
        c.insert(CobComp::new(set![3,4,5], set![12,13], 0, (0, 0)));

        assert_eq!(c.comps.len(), 2);
        
        c.insert(CobComp::new(set![2,3], set![11,12], 0, (0, 0)));

        assert_eq!(c.comps.len(), 1);
        assert_eq!(c, Cob::new(vec![
            CobComp::new(set![0,1,2,3,4,5], set![10,11,12,13], 0, (0, 0)),
        ]));
    }

    #[test]
    fn cob_comp_inv() { 
        let c0 = CobComp::new(set![0], set![1], 0, (0, 0));
        let c = Cob::new(vec![c0]);

        assert_eq!(c.is_invertible(), true);
        assert_eq!(c.inv(), Some(Cob::new(vec![
            CobComp::new(set![1], set![0], 0, (0, 0))
        ])));

        let c0 = CobComp::new(set![0], set![1], 0, (0, 0));
        let c1 = CobComp::new(set![2], set![3], 0, (0, 0));
        let c = Cob::new(vec![c0, c1]);

        assert_eq!(c.is_invertible(), true);
        assert_eq!(c.inv(), Some(Cob::new(vec![
            CobComp::new(set![1], set![0], 0, (0, 0)),
            CobComp::new(set![3], set![2], 0, (0, 0))
        ])));

        let c0 = CobComp::new(set![0], set![1], 0, (1, 0));
        let c = Cob::new(vec![c0]);

        assert_eq!(c.is_invertible(), false);
        assert_eq!(c.inv(), None);
    }
}