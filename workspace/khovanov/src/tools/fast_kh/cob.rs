use core::panic;
use std::hash::Hash;
use std::collections::HashSet;
use std::fmt::Display;
use derive_more::Display;
use itertools::Itertools;
use yui_core::{Elem, Ring, RingOps};
use yui_lin_comb::{FreeGen, OrdForDisplay};
use yui_link::Edge;
use yui_polynomial::Mono2;
use super::tng::{Tng, TngComp};

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Display)]
pub enum Dot { 
    None, X, Y
}

#[derive(PartialEq, Eq, Clone, Copy, Debug, Display)]
pub enum Bottom { 
    Src, Tgt
}

#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct CobComp { 
    src: Tng,
    tgt: Tng,
    genus: usize,
    dots: (usize, usize) // nums of X and Y dots resp. 
}

impl CobComp { 
    pub fn new(src: Tng, tgt: Tng, genus: usize, dots: (usize, usize)) -> Self { 
        debug_assert_eq!(src.endpts(), tgt.endpts());
        Self { src, tgt, genus, dots }
    }

    pub fn plain(src: Tng, tgt: Tng) -> Self { 
        Self::new(src, tgt, 0, (0, 0))
    }

    pub fn id(c: TngComp) -> Self { 
        Self::plain(
            Tng::from(c.clone()), 
            Tng::from(c),
        )
    }

    pub fn sdl(r0: (TngComp, TngComp), r1: (TngComp, TngComp)) -> Self { 
        assert!(r0.0.is_arc());
        assert!(r0.1.is_arc());
        assert!(r1.0.is_arc());
        assert!(r1.1.is_arc());
        assert!(r0.0 != r1.0);
        assert!(r0.0 != r1.1);
        assert!(r0.1 != r1.0);
        assert!(r0.1 != r1.1);

        Self::plain(
            Tng::new(vec![r0.0, r0.1]), 
            Tng::new(vec![r1.0, r1.1])
        )
    }

    pub fn cup(c: TngComp) -> Self { 
        assert!(c.is_circle());
        Self::plain(
            Tng::from(c),
            Tng::empty()
        )
    }

    pub fn cap(c: TngComp) -> Self { 
        Self::plain(
            Tng::empty(),
            Tng::from(c)
        )
    }

    pub fn genus(&self) -> usize { 
        self.genus
    }

    pub fn endpts(&self) -> HashSet<Edge> { 
        self.src.endpts() // == self.tgt.endpts()
    }

    pub fn bottom(&self, b: Bottom) -> &Tng { 
        match b { 
            Bottom::Src => &self.src,
            Bottom::Tgt => &self.tgt
        }
    }

    pub fn bottom_mut(&mut self, b: Bottom) -> &mut Tng { 
        match b { 
            Bottom::Src => &mut self.src,
            Bottom::Tgt => &mut self.tgt
        }
    }

    pub fn contains(&self, b: Bottom, c: &TngComp) -> bool { 
        self.bottom(b).contains(c)
    }

    pub fn index_of(&self, b: Bottom, c: &TngComp) -> Option<usize> { 
        self.bottom(b).index_of(c)
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

    pub fn is_cyl(&self) -> bool { // possibly with genus
        self.src.ncomps() == 1 && 
        self.tgt.ncomps() == 1
    }

    pub fn is_sdl(&self) -> bool { // possibly with genus
        self.src.ncomps() == 2 && 
        self.tgt.ncomps() == 2 && 
        self.src.comps().iter().all(|c| c.is_arc()) && 
        self.tgt.comps().iter().all(|c| c.is_arc()) && 
        self.src != self.tgt
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
        self.genus == 0 &&
        self.dots == (0, 0)
    }

    pub fn inv(&self) -> Option<Self> { 
        if self.is_invertible() { 
            let inv = Self::plain(
                self.tgt.clone(),
                self.src.clone() 
            );
            Some(inv)
        } else {
            None
        }
    }

    // χ(S) = 2 - 2g(S) - #(∂S)
    pub fn euler_num(&self) -> i32 { 
        let b = self.nbdr_comps() as i32;
        let g = self.genus as i32;
        2 - 2 * g - b
    }

    pub fn nbdr_comps(&self) -> usize { 
        let mut src_arcs: HashSet<_> = (0..self.src.ncomps()).filter(|&i| 
            self.src.comp(i).is_arc()
        ).collect();

        let mut tgt_arcs: HashSet<_> = (0..self.tgt.ncomps()).filter(|&i| 
            self.tgt.comp(i).is_arc()
        ).collect();

        assert_eq!(src_arcs.len(), tgt_arcs.len());

        let src_circs = self.src.ncomps() - src_arcs.len();
        let tgt_circs = self.tgt.ncomps() - tgt_arcs.len();

        let mut side_circs = 0;

        while !src_arcs.is_empty() { 
            let mut i0 = src_arcs.iter().next().cloned().unwrap();
            loop { 
                src_arcs.remove(&i0);

                let c0 = self.src.comp(i0);
                let Some(j) = tgt_arcs.iter().find(|&&j| { 
                    self.tgt.comp(j).is_connectable(&c0)
                }).cloned() else { panic!() };

                tgt_arcs.remove(&j);

                let c1 = self.tgt.comp(j);
                if let Some(i1) = src_arcs.iter().find(|&&i| { 
                    i != i0 && self.src.comp(i).is_connectable(&c1)
                }).cloned() { 
                    i0 = i1;
                } else { 
                    side_circs += 1;
                    break
                }
            }
        }

        src_circs + tgt_circs + side_circs
    }

    pub fn add_dot(&mut self, dot: Dot) { 
        match dot { 
            Dot::X => self.dots.0 += 1,
            Dot::Y => self.dots.1 += 1,
            _      => ()
        }
    }

    pub fn cap_off(&mut self, b: Bottom, i: usize) {
        assert!(self.bottom(b).comp(i).is_circle());
        self.bottom_mut(b).remove_at(i);
    }

    // connect = horizontal composition
    pub fn is_connectable(&self, other: &Self) -> bool { 
        self.src.comps().iter().any(|c1| 
            if c1.is_arc() { 
                other.src.comps().iter().any(|c2| { 
                    c2.is_arc() && c1.is_connectable(c2)
                })
            } else { 
                false 
            }
        )
    }
    
    pub fn connect(&mut self, other: Self) { 
        debug_assert!(self.is_connectable(&other));

        // χ(S∪S') = χ(S) + χ(S') - χ(S∩S'),
        // χ(S) = 2 - 2g(S) - #∂S, 
        // χ(S∩S') = #(S∩S').
        // → g(S∪S') = g(S) + g(S') + 1/2 #{∂S + ∂S' + S∩S' - ∂(S∪S')} - 1.

        let g1 = self.genus;
        let g2 = other.genus;
        let b1 = self.nbdr_comps();
        let b2 = other.nbdr_comps();
        let f = self.endpts().intersection(&other.endpts()).count();

        let CobComp{ src, tgt, genus: _, dots } = other;

        self.src.connect(src);
        self.tgt.connect(tgt);

        let b = self.nbdr_comps();
        assert!((b1 + b2 + f - b) % 2 == 0);
        
        self.genus = g1 + g2 + (b1 + b2 + f - b)/2 - 1;

        self.dots.0 += dots.0;
        self.dots.1 += dots.1;
    }

    pub fn eval<R>(&self, h: &R, t: &R) -> R
    where R: Ring, for<'x> &'x R: RingOps<R> {
        assert!(self.is_closed());

        // TODO must `neck cut` until g = 0!

        fn eval<R>(x: usize, y: usize, h: &R, t: &R) -> R
        where R: Ring, for<'x> &'x R: RingOps<R> { 
            match (x, y) { 
                (0, 0) => R::zero(),
                (1, 0) | (0, 1)  => R::one(),
                (x, y) if x >= 1 && y >= 1 => // XY = t
                    t * eval(x-1, y-1, h, t),
                (x, 0) if x >= 2 => // X^2 = hX + t
                    h * eval(x-1, 0, h, t) + 
                    t * eval(x-2, 0, h, t),
                (0, y) if y >= 2 => // Y^2 = -hY + t
                    -h * eval(0, y-1, h, t) + 
                     t * eval(0, y-2, h, t),
                _ => panic!()
            }
        }

        let (x, y) = self.dots;
        eval(x, y, h, t)
    }
}

impl Display for CobComp {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match (self.src.ncomps(), self.tgt.ncomps()) { 
            (0, 0) => "S",
            (0, 1) => "ι",
            (1, 0) => "ε",
            (1, 1) => "cyl",
            (2, 1) => "m",
            (1, 2) => "Δ",
            (2, 2) if self.is_sdl() => "sdl",
            _      => "cob"
        };

        let g = if self.genus > 0 { 
            yui_utils::subscript(self.genus as isize)
        } else { 
            "".to_string()
        };

        let dots = if self.has_dots() { 
            let (p, q) = self.dots;
            format!("({})", Mono2::<'X','Y', _>::new(p, q).to_string())
        } else { 
            String::new()
        };
        
        let cob = format!("{s}{g}{dots}");

        if self.src.is_empty() && self.tgt.is_empty() { 
            write!(f, "{cob}")
        } else { 
            write!(f, "{cob}({} -> {})", &self.src, &self.tgt)
        }
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
    
    pub fn id(v: &Tng) -> Self { 
        let comps = (0..v.ncomps()).map(|i| {
            let c = v.comp(i).clone();
            CobComp::id(c)
        }).collect();
        Self::new(comps)
    }

    pub fn sdl(src: &Tng, tgt: &Tng) -> Self { 
        assert_eq!(src.ncomps(), 2);
        assert_eq!(tgt.ncomps(), 2);

        let r0 = (src.comp(0).clone(), src.comp(1).clone());
        let r1 = (tgt.comp(0).clone(), tgt.comp(1).clone());
        let sdl = CobComp::sdl(r0, r1);

        Self::from(sdl)
    }

    pub fn ncomps(&self) -> usize { 
        self.comps.len()
    }

    pub fn comp(&self, i: usize) -> &CobComp { 
        &self.comps[i]
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

    pub fn euler_num(&self) -> i32 { 
        self.comps.iter().map(|c| c.euler_num()).sum()
    }

    pub fn nbdr_comps(&self) -> usize { 
        self.comps.iter().map(|c| c.nbdr_comps()).sum()
    }

    pub fn cap_off(&mut self, b: Bottom, c: &TngComp, x: Dot) {
        assert!(c.is_circle());
        let Some((i, comp, p)) = self.find_comp(b, c) else { 
            panic!("{c} not found in {} ({b})", self)
        };

        comp.cap_off(b, p);
        comp.add_dot(x);

        if comp.is_removable() { 
            self.comps.remove(i);
        }
    }

    fn find_comp(&mut self, b: Bottom, c: &TngComp) -> Option<(usize, &mut CobComp, usize)> { 
        self.comps.iter_mut().enumerate().filter_map(|(i, comp)| 
            if let Some(p) = comp.index_of(b, c) { 
                Some((i, comp, p))
            } else { 
                None
            }
        ).next()
    }

    pub fn connect(&mut self, cob: Cob) { // horizontal composition
        for c in cob.comps.into_iter() { 
            self.connect_comp(c);
        }
    }

    pub fn connect_comp(&mut self, mut c: CobComp) {
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
        } else if self.is_zero() { 
            write!(f, "0")
        } else if self.comps.len() == 1 {
            self.comps[0].fmt(f)
        } else { 
            let cobs = self.comps.iter().map(|c| c.to_string()).join(" ⊔ ");
            write!(f, "{{{}}}", cobs)
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
    fn cob_contains() { 
        let src = Tng::new(vec![
            TngComp::arc(1,2),
            TngComp::arc(3,4),
            TngComp::circ(5),
        ]);
        let tgt = Tng::new(vec![
            TngComp::arc(1,3),
            TngComp::arc(2,4),
            TngComp::circ(6),
        ]);
        let c = CobComp::plain(src, tgt);
        
        let c0 = TngComp::arc(1,2);
        let c1 = TngComp::circ(6);

        assert!( c.contains(Bottom::Src, &c0));
        assert!(!c.contains(Bottom::Src, &c1));
        assert!(!c.contains(Bottom::Tgt, &c0));
        assert!( c.contains(Bottom::Tgt, &c1));
    }

    #[test]
    fn is_connectable() { 
        let src = Tng::new(vec![
            TngComp::arc(1,2),
            TngComp::arc(3,4),
            TngComp::circ(10),
        ]);
        let tgt = Tng::new(vec![
            TngComp::arc(1,3),
            TngComp::arc(2,4),
            TngComp::circ(11),
        ]);
        let c = CobComp::plain(src, tgt);

        let c1 = CobComp::id(
            TngComp::arc(0,1)
        );
        let c2 = CobComp::sdl(
            (TngComp::arc(0,1), TngComp::arc(90,91)),
            (TngComp::arc(0,90), TngComp::arc(1,91)),
        );
        let c3 = CobComp::id(TngComp::arc(5,6));

        assert!(c.is_connectable(&c1));
        assert!(c.is_connectable(&c2));
        assert!(!c.is_connectable(&c3));
    }

    #[test]
    fn connect1() { 
        let src = Tng::new(vec![
            TngComp::arc(1,2),
            TngComp::arc(3,4),
            TngComp::circ(10),
        ]);
        let tgt = Tng::new(vec![
            TngComp::arc(1,3),
            TngComp::arc(2,4),
            TngComp::circ(11),
        ]);

        let mut c = CobComp::plain(src, tgt);
        c.connect(CobComp::id(
            TngComp::arc(0,1)
        ));

        assert_eq!(c, CobComp::plain(
            Tng::new(vec![
                TngComp::arc(0,2), // [0,1,2] -> [0,2]
                TngComp::arc(3,4),
                TngComp::circ(10),
            ]),
            Tng::new(vec![
                TngComp::arc(0,3), // [0,1,2] -> [0,2]
                TngComp::arc(2,4),
                TngComp::circ(11),
            ])
        ));
    }

    #[test]
    fn connect2() { 
        let src = Tng::new(vec![
            TngComp::arc(1,2),
            TngComp::arc(3,4),
            TngComp::circ(10),
        ]);
        let tgt = Tng::new(vec![
            TngComp::arc(1,3),
            TngComp::arc(2,4),
            TngComp::circ(11),
        ]);

        let mut c = CobComp::plain(src, tgt);
        c.connect(CobComp::id(
            TngComp::arc(1,3)
        ));

        assert_eq!(c, CobComp::plain(
            Tng::new(vec![
                TngComp::arc(4,2),
                TngComp::circ(10),
            ]),
            Tng::new(vec![
                TngComp::circ(1),
                TngComp::arc(2,4),
                TngComp::circ(11),
            ])
        ));
    }

    #[test]
    fn euler_num() { 
        let c0 = CobComp::id(
            TngComp::arc(1,2)
        );
        let c1 = CobComp::sdl(
            (TngComp::arc(3,4), TngComp::arc(5,6)),
            (TngComp::arc(4,5), TngComp::arc(6,3)),
        );
        let c2 = CobComp::plain(
            Tng::from(TngComp::circ(10)),
            Tng::new(vec![TngComp::circ(10), TngComp::circ(11)]),
        );
        let c3 = CobComp::cup(
            TngComp::circ(20)
        );
        let c4 = CobComp::cap(
            TngComp::circ(30)
        );

        assert_eq!(c0.nbdr_comps(), 1);
        assert_eq!(c1.nbdr_comps(), 1);
        assert_eq!(c2.nbdr_comps(), 3);
        assert_eq!(c3.nbdr_comps(), 1);
        assert_eq!(c4.nbdr_comps(), 1);

        assert_eq!(c0.euler_num(), 1);
        assert_eq!(c1.euler_num(), 1);
        assert_eq!(c2.euler_num(), -1);
        assert_eq!(c3.euler_num(), 1);
        assert_eq!(c4.euler_num(), 1);

        let cob = Cob::new(vec![c0,c1,c2,c3,c4]);
        assert_eq!(cob.euler_num(), 3);
        assert_eq!(cob.nbdr_comps(), 7);
    }

    #[test]
    fn connect_incr_genus() { 
        let mut c0 = CobComp::plain(
            Tng::new(vec![
                TngComp::arc(1,2),
                TngComp::arc(3,4)
            ]),
            Tng::new(vec![
                TngComp::arc(1,2),
                TngComp::arc(3,4)
            ]),
        );
        let c1 = CobComp::id(
            TngComp::arc(1, 3)
        );
        let c2 = CobComp::id(
            TngComp::arc(2, 4)
        );

        assert_eq!(c0.genus, 0);
        assert_eq!(c1.genus, 0);
        assert_eq!(c2.genus, 0);

        c0.connect(c1);

        assert_eq!(c0.genus, 1);
        assert_eq!(c0.euler_num(), -1);

        c0.connect(c2);

        assert_eq!(c0.genus, 1);
        assert_eq!(c0.euler_num(), -2);

        c0.cap_off(Bottom::Src, 0);

        assert_eq!(c0.genus, 1);
        assert_eq!(c0.euler_num(), -1);

        c0.cap_off(Bottom::Tgt, 0);

        assert_eq!(c0.genus, 1);
        assert_eq!(c0.euler_num(), 0);
        assert!(c0.is_closed()); // torus
    }

    #[test]
    fn inv() { 
        let cc0 = CobComp::id(TngComp::arc(0, 1));
        let cc1 = CobComp::plain(
            Tng::from(TngComp::circ(2)), 
            Tng::from(TngComp::circ(3))
        );

        assert_eq!(cc0.is_invertible(), true);
        assert_eq!(cc0.inv(), Some(cc0.clone()));

        assert_eq!(cc1.is_invertible(), true);
        assert_eq!(cc1.inv(), Some(CobComp::plain(
            Tng::from(TngComp::circ(3)), 
            Tng::from(TngComp::circ(2))
        )));

        let c0 = Cob::new(vec![cc0, cc1]);

        assert_eq!(c0.is_invertible(), true);
        assert_eq!(c0.inv(), Some(Cob::new(vec![
            c0.comp(0).inv().unwrap(),
            c0.comp(1).inv().unwrap()
        ])));

        let c1 = Cob::from(
            CobComp::sdl(
                (TngComp::arc(1,2), TngComp::arc(3,4)),
                (TngComp::arc(1,3), TngComp::arc(2,4))
            )
        );

        assert_eq!(c1.is_invertible(), false);
        assert_eq!(c1.inv(), None);

        let mut c2 = c0.clone();
        c2.comps[0].add_dot(Dot::X);

        assert_eq!(c2.is_invertible(), false);
        assert_eq!(c2.inv(), None);

        let mut c3 = c0.clone();
        c3.comps[0].genus += 1;

        assert_eq!(c3.is_invertible(), false);
        assert_eq!(c3.inv(), None);
    }
}