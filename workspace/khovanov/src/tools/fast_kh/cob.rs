use core::panic;
use std::hash::Hash;
use std::collections::HashSet;
use std::fmt::Display;
use derive_more::Display;
use itertools::Itertools;
use yui_core::{Elem, Ring, RingOps};
use yui_lin_comb::{FreeGen, OrdForDisplay};
use yui_link::{Component, Edge};
use yui_polynomial::Mono2;
use super::tng::Tng;

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Display)]
pub enum Dot { 
    None, X, Y
}

#[derive(PartialEq, Eq, Clone, Copy, Debug)]
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

    pub fn new_plain(src: Tng, tgt: Tng) -> Self { 
        Self::new(src, tgt, 0, (0, 0))
    }

    pub fn id(c: Component) -> Self { 
        assert!(!c.is_empty());
        Self::new_plain(
            Tng::from(c.clone()), 
            Tng::from(c),
        )
    }

    pub fn sdl(r0: (Component, Component), r1: (Component, Component)) -> Self { 
        // TODO validate
        Self::new_plain(
            Tng::new(vec![r0.0, r0.1]), 
            Tng::new(vec![r1.0, r1.1])
        )
    }

    pub fn cup(c: Component) -> Self { 
        assert!(c.is_circle());
        Self::new_plain(
            Tng::from(c),
            Tng::empty()
        )
    }

    pub fn cap(c: Component) -> Self { 
        Self::new_plain(
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

    pub fn contains(&self, c: &Component, e: End) -> bool { 
        self.end(e).contains(c)
    }

    pub fn index_of(&self, c: &Component, e: End) -> Option<usize> { 
        self.end(e).index_of(c)
    }

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
        self.src.ncomps() == 1 && 
        self.tgt.ncomps() == 1 && 
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
            let inv = Self::new_plain(
                self.tgt.clone(),
                self.src.clone() 
            );
            Some(inv)
        } else {
            None
        }
    }

    pub fn cap_off(&mut self, i: usize, e: End) {
        assert!(self.end(e).comp(i).is_circle());
        self.end_mut(e).remove_at(i);
    }

    pub fn add_dot(&mut self, dot: Dot) { 
        match dot { 
            Dot::X => self.dots.0 += 1,
            Dot::Y => self.dots.1 += 1,
            _      => ()
        }
    }

    fn end(&self, e: End) -> &Tng { 
        if e.is_src() { 
            &self.src 
        } else { 
            &self.tgt
        }
    }

    fn end_mut(&mut self, e: End) -> &mut Tng { 
        if e.is_src() { 
            &mut self.src 
        } else { 
            &mut self.tgt
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
        
        write!(f, "{s}({} -> {}{})", &self.src, &self.tgt, dots)
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
    
    pub fn id(v: &Tng) -> Self { 
        let comps = (0..v.ncomps()).map(|i| {
            let c = v.comp(i).clone();
            CobComp::id(c)
        }).collect();
        Self::new(comps)
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

    pub fn cap_off(&mut self, c: &Component, x: Dot, e: End) {
        assert!(c.is_circle());
        let Some((i, comp, p)) = self.find_comp(c, e) else { 
            panic!("{c} not found in {} ({e:?})", self)
        };

        comp.cap_off(p, e);
        comp.add_dot(x);

        if comp.is_removable() { 
            self.comps.remove(i);
        }
    }

    fn find_comp(&mut self, c: &Component, e: End) -> Option<(usize, &mut CobComp, usize)> { 
        self.comps.iter_mut().enumerate().filter_map(|(i, comp)| 
            if let Some(p) = comp.index_of(c, e) { 
                Some((i, comp, p))
            } else { 
                None
            }
        ).next()
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
    use yui_link::Component;

    use super::CobComp;
    use super::*;
 
    #[test]
    fn cob_contains() { 
        let src = Tng::new(vec![
            Component::arc(vec![1,2]),
            Component::arc(vec![3,4]),
            Component::circ(vec![5]),
        ]);
        let tgt = Tng::new(vec![
            Component::arc(vec![1,3]),
            Component::arc(vec![2,4]),
            Component::circ(vec![6]),
        ]);
        let c = CobComp::new_plain(src, tgt);
        
        let c0 = Component::arc(vec![1,2]);
        let c1 = Component::circ(vec![6]);

        assert!( c.contains(&c0, End::Src));
        assert!(!c.contains(&c1, End::Src));
        assert!(!c.contains(&c0, End::Tgt));
        assert!( c.contains(&c1, End::Tgt));
    }

    #[test]
    fn is_connectable() { 
        let src = Tng::new(vec![
            Component::arc(vec![1,2]),
            Component::arc(vec![3,4]),
            Component::circ(vec![10]),
        ]);
        let tgt = Tng::new(vec![
            Component::arc(vec![1,3]),
            Component::arc(vec![2,4]),
            Component::circ(vec![11]),
        ]);
        let c = CobComp::new_plain(src, tgt);

        let c1 = CobComp::id(
            Component::arc(vec![0,1])
        );
        let c2 = CobComp::sdl(
            (Component::arc(vec![0,1]), Component::arc(vec![90,91])),
            (Component::arc(vec![0,90]), Component::arc(vec![1,91])),
        );
        let c3 = CobComp::id(Component::arc(vec![5,6]));

        assert!(c.is_connectable(&c1));
        assert!(c.is_connectable(&c2));
        assert!(!c.is_connectable(&c3));
    }

    #[test]
    fn connect1() { 
        let src = Tng::new(vec![
            Component::arc(vec![1,2]),
            Component::arc(vec![3,4]),
            Component::circ(vec![10]),
        ]);
        let tgt = Tng::new(vec![
            Component::arc(vec![1,3]),
            Component::arc(vec![2,4]),
            Component::circ(vec![11]),
        ]);

        let mut c = CobComp::new_plain(src, tgt);
        c.connect(CobComp::id(
            Component::arc(vec![0,1])
        ));

        assert_eq!(c, CobComp::new_plain(
            Tng::new(vec![
                Component::arc(vec![0,1,2]),
                Component::arc(vec![3,4]),
                Component::circ(vec![10]),
            ]),
            Tng::new(vec![
                Component::arc(vec![0,1,3]),
                Component::arc(vec![2,4]),
                Component::circ(vec![11]),
            ])
        ));
    }

    #[test]
    fn connect2() { 
        let src = Tng::new(vec![
            Component::arc(vec![1,2]),
            Component::arc(vec![3,4]),
            Component::circ(vec![10]),
        ]);
        let tgt = Tng::new(vec![
            Component::arc(vec![1,3]),
            Component::arc(vec![2,4]),
            Component::circ(vec![11]),
        ]);

        let mut c = CobComp::new_plain(src, tgt);
        c.connect(CobComp::id(
            Component::arc(vec![1,3])
        ));

        assert_eq!(c, CobComp::new_plain(
            Tng::new(vec![
                Component::arc(vec![4,3,1,2]),
                Component::circ(vec![10]),
            ]),
            Tng::new(vec![
                Component::circ(vec![1,3]),
                Component::arc(vec![2,4]),
                Component::circ(vec![11]),
            ])
        ));
    }

    #[test]
    fn euler_num() { 
        let c0 = CobComp::id(
            Component::arc(vec![1,2])
        );
        let c1 = CobComp::sdl(
            (Component::arc(vec![3,4]), Component::arc(vec![5,6])),
            (Component::arc(vec![4,5]), Component::arc(vec![6,3])),
        );
        let c2 = CobComp::new_plain(
            Tng::from(Component::circ(vec![10])),
            Tng::new(vec![Component::circ(vec![10]), Component::circ(vec![11])]),
        );
        let c3 = CobComp::cup(
            Component::circ(vec![20])
        );
        let c4 = CobComp::cap(
            Component::circ(vec![30])
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
        let mut c0 = CobComp::new_plain(
            Tng::new(vec![
                Component::arc(vec![1,2]),
                Component::arc(vec![3,4])
            ]),
            Tng::new(vec![
                Component::arc(vec![1,2]),
                Component::arc(vec![3,4])
            ]),
        );
        let c1 = CobComp::id(
            Component::arc(vec![1, 3])
        );
        let c2 = CobComp::id(
            Component::arc(vec![2, 4])
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

        c0.cap_off(0, End::Src);

        assert_eq!(c0.genus, 1);
        assert_eq!(c0.euler_num(), -1);

        c0.cap_off(0, End::Tgt);

        assert_eq!(c0.genus, 1);
        assert_eq!(c0.euler_num(), 0);
        assert!(c0.is_closed()); // torus
    }


    // #[test]
    // fn cob_comp_inv() { 
    //     let c0 = CobComp::new(set![0], set![1], 0, (0, 0));
    //     let c = Cob::new(vec![c0]);

    //     assert_eq!(c.is_invertible(), true);
    //     assert_eq!(c.inv(), Some(Cob::new(vec![
    //         CobComp::new(set![1], set![0], 0, (0, 0))
    //     ])));

    //     let c0 = CobComp::new(set![0], set![1], 0, (0, 0));
    //     let c1 = CobComp::new(set![2], set![3], 0, (0, 0));
    //     let c = Cob::new(vec![c0, c1]);

    //     assert_eq!(c.is_invertible(), true);
    //     assert_eq!(c.inv(), Some(Cob::new(vec![
    //         CobComp::new(set![1], set![0], 0, (0, 0)),
    //         CobComp::new(set![3], set![2], 0, (0, 0))
    //     ])));

    //     let c0 = CobComp::new(set![0], set![1], 0, (1, 0));
    //     let c = Cob::new(vec![c0]);

    //     assert_eq!(c.is_invertible(), false);
    //     assert_eq!(c.inv(), None);
    // }
}