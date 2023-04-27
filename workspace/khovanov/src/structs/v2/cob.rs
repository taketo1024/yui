use core::panic;
use std::hash::Hash;
use std::collections::{HashSet, VecDeque};
use std::fmt::Display;
use std::ops::{Mul, MulAssign};
use auto_impl_ops::auto_ops;
use derive_more::Display;
use itertools::Itertools;
use num_traits::Zero;
use cartesian::cartesian;
use yui_core::{Elem, Ring, RingOps};
use yui_lin_comb::{FreeGen, OrdForDisplay, LinComb};
use yui_link::{Edge, Crossing, Resolution};
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

    pub fn merge(from: (TngComp, TngComp), to: TngComp) -> Self { 
        assert!(from.0.is_circle() || from.1.is_circle());
        Self::plain(
            Tng::new(vec![from.0, from.1]), 
            Tng::new(vec![to])
        )
    }

    pub fn split(from: TngComp, to: (TngComp, TngComp)) -> Self { 
        assert!(to.0.is_circle() || to.1.is_circle());
        Self::plain(
            Tng::new(vec![from]), 
            Tng::new(vec![to.0, to.1])
        )
    }

    pub fn cup(c: TngComp) -> Self { 
        assert!(c.is_circle());
        Self::plain(
            Tng::empty(),
            Tng::from(c)
        )
    }

    pub fn cap(c: TngComp) -> Self { 
        Self::plain(
            Tng::from(c),
            Tng::empty()
        )
    }

    pub fn closed(g: usize) -> Self { 
        Self::new(
            Tng::empty(),
            Tng::empty(),
            g,
            (0, 0)
        )
    }

    pub fn sphere() -> Self { 
        Self::closed(0)
    }

    pub fn src(&self) -> &Tng { 
        &self.src
    }

    pub fn tgt(&self) -> &Tng { 
        &self.tgt
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

    pub fn is_id(&self) -> bool { 
        self.is_cyl() && 
        self.src.comp(0) == self.tgt.comp(0) && 
        self.genus == 0
    }

    pub fn is_sdl(&self) -> bool { // possibly with genus
        self.src.ncomps() == 2 && 
        self.tgt.ncomps() == 2 && 
        self.src.comps().iter().all(|c| c.is_arc()) && 
        self.tgt.comps().iter().all(|c| c.is_arc()) && 
        self.src != self.tgt
    }

    pub fn is_zero(&self) -> bool { 
        self.is_closed() && 
        self.genus % 2 == 0 &&
        self.dots.0 == self.dots.1 // XY = T
    }

    pub fn is_removable(&self) -> bool { // TODO rename to `is_one()`
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

        // χ(S∪S') = χ(S) + χ(S') - χ(S∩S')
        //         = 2 - 2g(S∪S') - #∂(S∪S'),
        // χ(S∩S') = #{ arcs in S∩S' }.
        // 2g(S∪S') = 2 - (χ(S) + χ(S') + #∂(S∪S')) + #∂(S∩S').

        let x1 = self.euler_num();
        let x2 = other.euler_num();

        let a = self.endpts().intersection(
            &other.endpts()
        ).count() as i32;

        assert!(a > 0);

        let CobComp{ src, tgt, genus: _, dots } = other;

        self.src.connect(src);
        self.tgt.connect(tgt);

        let b = self.nbdr_comps() as i32;
        let g = 2 - (x1 + x2 + b) + a;

        assert!(g >= 0);
        assert!(g % 2 == 0);
        
        self.genus = (g / 2) as usize;

        self.dots.0 += dots.0;
        self.dots.1 += dots.1;
    }

    pub fn should_part_eval(&self) -> bool {
        self.is_zero() ||
        self.genus > 0 || 
        self.dots.0 >= 1 && self.dots.1 >= 1 ||
        self.dots.0 >= 2 ||
        self.dots.1 >= 2
    }

    pub fn part_eval<R>(&self, h: &R, t: &R) -> LinComb<CobComp, R>
    where R: Ring, for<'x> &'x R: RingOps<R> {

        fn cob(c: &CobComp, x: usize, y: usize) -> CobComp { 
            CobComp { src: c.src.clone(), tgt: c.tgt.clone(), genus: 0, dots: (x, y) }
        }

        fn eval<R>(c: &CobComp, g: usize, x: usize, y: usize, h: &R, t: &R) -> LinComb<CobComp, R>
        where R: Ring, for<'x> &'x R: RingOps<R> { 
            match (g, x, y) { 
                (g, _, _) if g > 0 => { // neck-cut
                    eval(c, g-1, x+1, y, h, t) + 
                    eval(c, g-1, x, y+1, h, t)
                },
                (0, x, y) if x >= 1 && y >= 1 => // XY = t
                    eval(c, 0, x-1, y-1, h, t) * t,
                (0, x, 0) if x >= 2 => // X^2 = hX + t
                    eval(c, 0, x-1, 0, h, t) * h + 
                    eval(c, 0, x-2, 0, h, t) * t,
                (0, 0, y) if y >= 2 => // Y^2 = -hY + t
                    eval(c, 0, 0, y-1, h, t) * -h + 
                    eval(c, 0, 0, y-2, h, t) *  t,
                (0, 0, 0) if c.is_closed() => // ε.ι = 0
                    LinComb::zero(),
                _ =>
                    LinComb::from_gen(cob(c, x, y))
            }
        }

        let g = self.genus;
        let (x, y) = self.dots;

        eval(self, g, x, y, h, t)
    }

    pub fn eval<R>(&self, h: &R, t: &R) -> R
    where R: Ring, for<'x> &'x R: RingOps<R> {
        assert!(self.is_closed());

        self.part_eval(h, t).into_iter().map(|(c, r)| {
            debug_assert!(c.is_removable()); // either ε.X.ι or ε.Y.ι .
            r
        }).sum()
    }
}

impl From<&Crossing> for CobComp {
    fn from(x: &Crossing) -> Self {
        let src = Tng::res(x, Resolution::Res0);
        let tgt = Tng::res(x, Resolution::Res1);
        Self::plain(src, tgt)
    }
}

impl Display for CobComp {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match (self.src.ncomps(), self.tgt.ncomps()) { 
            (0, 0) => "S",
            (0, 1) => "ι",
            (1, 0) => "ε",
            (1, 1) if self.is_id() => "id",
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
            format!("{}.", Mono2::<'X','Y', _>::new(p, q).to_string())
        } else { 
            String::new()
        };
        
        let cob = format!("{dots}{s}{g}");

        if self.src.is_empty() && self.tgt.is_empty() { 
            write!(f, "{cob}")
        } else if self.is_id() { 
            write!(f, "{cob}({})", &self.src.comp(0))
        } else { 
            write!(f, "{cob}({} -> {})", &self.src, &self.tgt)
        }
    }
}

impl PartialOrd for CobComp {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for CobComp {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        use std::cmp::Ordering::*;

        fn cmp_opt(c0: Option<&TngComp>, c1: Option<&TngComp>) -> Option<std::cmp::Ordering> {
            match (c0, c1) { 
                (Some(c0), Some(c1)) => Some(c0.cmp(c1)),
                (None, Some(_)) => Some(Less),
                (Some(_), None) => Some(Greater),
                (None, None) => None
            }
        }

        if let Some(src_cmp) = cmp_opt(
            self.src.min_comp(), 
            other.src.min_comp()
        ) { 
            src_cmp
        } else if let Some(tgt_cmp) = cmp_opt(
            self.tgt.min_comp(), 
            other.tgt.min_comp()
        ) { 
            tgt_cmp
        } else { 
            [self.genus, self.dots.0, self.dots.1].cmp(
                &[other.genus, other.dots.0, other.dots.1]
            )
        }
    }
}

impl Default for CobComp {
    fn default() -> Self {
        Self::sphere() // == zero
    }
}

impl Elem for CobComp {
    fn set_symbol() -> String {
        "CobComp".to_string()
    }
}

impl FreeGen for CobComp {}

#[derive(Clone, PartialEq, Eq, Hash, Debug, Default)]
pub struct Cob { 
    comps: Vec<CobComp>
}

impl Cob {
    pub fn new(comps: Vec<CobComp>) -> Self { 
        let mut cob = Self { comps };
        cob.sort_comps();
        cob
    }

    pub fn empty() -> Self { 
        Self::new(vec![])
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

    pub fn src(&self) -> Tng { 
        self.comps.iter().fold(Tng::empty(), |mut t, c| {
            t.connect(c.src.clone());
            t
        })
    }

    pub fn tgt(&self) -> Tng { 
        self.comps.iter().fold(Tng::empty(), |mut t, c| {
            t.connect(c.tgt.clone());
            t
        })
    }

    pub fn is_empty(&self) -> bool { 
        self.comps.is_empty()
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

        self.sort_comps()
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

    pub fn connect(&mut self, other: Cob) { // horizontal composition
        for c in other.comps.into_iter() { 
            self.connect_comp(c);
        }
        self.sort_comps()
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

    pub fn is_stackable(&self, other: &Self) -> bool { 
        self.comps.iter().fold(0, |n, c| n + c.tgt.ncomps()) == 
        other.comps.iter().fold(0, |n, c| n + c.src.ncomps()) && 
        self.comps.iter().all(|c| c.tgt.comps().iter().all(|a|
            other.comps.iter().any(|c| c.contains(Bottom::Src, &a))
        ))
    }

    pub fn stack(&mut self, other: Cob) { // vertical composition
        debug_assert!(self.is_stackable(&other));

        if self.is_empty() { 
            *self = other;
            return;
        } else if other.is_empty() { 
            return;
        }

        let mut bot = std::mem::replace(&mut self.comps, vec![]);
        let mut top = other.comps;

        while !(bot.is_empty() && top.is_empty()) { 
            let (mut b, mut t) = self.take_stackable_comps(&mut bot, &mut top);

            if t.is_empty() { 
                assert_eq!(b.len(), 1);
                self.comps.push(b.remove(0))
            } else if b.is_empty() {
                assert_eq!(t.len(), 1);
                self.comps.push(t.remove(0))
            } else { 
                let c = self.stack_comps(b, t);
                self.comps.push(c)
            }
        }

        self.sort_comps()
    }

    // collect `bot` & `top` comps that will form a connected component.
    fn take_stackable_comps(&self, bot: &mut Vec<CobComp>, top: &mut Vec<CobComp>) -> (Vec<CobComp>, Vec<CobComp>) {
        let mut res_bot = vec![];
        let mut res_top = vec![];

        let mut q_bot = VecDeque::new();
        let mut q_top = VecDeque::new();

        if !bot.is_empty() { 
            q_bot.push_back(bot.remove(0))
        } else if !top.is_empty() { 
            q_top.push_back(top.remove(0))
        }

        while !(q_bot.is_empty() && q_top.is_empty()) {
            while let Some(b) = q_bot.pop_front() {
                for c in b.tgt.comps().iter() { 
                    if let Some(i) = top.iter().position(|t| t.src.contains(&c)) {
                        let t = top.remove(i);
                        q_top.push_back(t);
                    }
                }
                res_bot.push(b)
            }
            while let Some(t) = q_top.pop_front() {
                for c in t.src.comps().iter() { 
                    if let Some(i) = bot.iter().position(|b| b.tgt.contains(&c)) {
                        let b = bot.remove(i);
                        q_bot.push_back(b);
                    }
                }
                res_top.push(t)
            }
        }

        (res_bot, res_top)
    }

    fn stack_comps(&self, bot: Vec<CobComp>, top: Vec<CobComp>) -> CobComp { 
        // `bot`, `top` must be non-empty for genus calculation. 
        assert!(!bot.is_empty());
        assert!(!top.is_empty());

        let x0: i32 = bot.iter().map(|c| c.euler_num()).sum();
        let x1: i32 = top.iter().map(|c| c.euler_num()).sum();
        
        let a : i32 = bot.iter().map(|c| 
            c.tgt.comps().iter().filter(|a| a.is_arc()).count() as i32
        ).sum();

        let dots = bot.iter().chain(top.iter()).fold((0, 0), |mut res, c| { 
            res.0 += c.dots.0;
            res.1 += c.dots.1;
            res
        });

        let src = bot.into_iter().fold(Tng::empty(), |mut res, c| {
            res.connect(c.src);
            res
        });

        let tgt = top.into_iter().fold(Tng::empty(), |mut res, c| {
            res.connect(c.tgt);
            res
        });

        let mut c = CobComp::new(src, tgt, 0, dots);
        let b = c.nbdr_comps() as i32;
        let g = 2 - (x0 + x1 + b) + a;

        assert!(g >= 0);
        assert!(g % 2 == 0);
        
        c.genus = (g / 2) as usize;

        c
    }

    pub fn should_part_eval(&self) -> bool {
        self.comps.iter().any(|c| c.should_part_eval())
    }

    pub fn part_eval<R>(self, h: &R, t: &R) -> LinComb<Cob, R>
    where R: Ring, for<'x> &'x R: RingOps<R> {
        if self.should_part_eval() { 
            let init = LinComb::from_gen(Cob::empty());
            self.comps.iter().fold(init, |res, c| { 
                let eval = c.part_eval(h, t);
                let all = cartesian!(res.iter(), eval.iter());
                let prod = all.map(|((cob, r), (c, s))| {
                    let mut cob = cob.clone();
                    cob.comps.push(c.clone());
                    (cob, r * s)
                }).map(|(mut cob, r)| {
                    cob.sort_comps();
                    (cob, r)
                });
                prod.collect()
            })
        } else { 
            LinComb::from_gen(self)
        }
    }

    pub fn eval<R>(&self, h: &R, t: &R) -> R
    where R: Ring, for<'x> &'x R: RingOps<R> {
        self.comps.iter().map(|c| 
            c.eval(h, t)
        ).product()
    }

    fn sort_comps(&mut self) { 
        self.comps.sort()
    }

    fn min_comp(&self) -> Option<&CobComp> { 
        self.comps.iter().min_by(|c0, c1| c0.cmp_for_display(c1))
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
            write!(f, "|{}|", cobs)
        }
    }
}

impl OrdForDisplay for Cob {
    fn cmp_for_display(&self, other: &Self) -> std::cmp::Ordering {
        use std::cmp::Ordering::*;

        let Some(c0) = self.min_comp() else {
            return Less;
        };
        let Some(c1) = other.min_comp() else {
            return Greater;
        };
        c0.cmp(c1)
    }
}

impl Elem for Cob {
    fn set_symbol() -> String {
        "Cob".to_string()
    }
}

impl FreeGen for Cob {}

#[auto_ops]
impl Mul for Cob {
    type Output = Cob;
    fn mul(self, mut rhs: Self) -> Self::Output {
        rhs.stack(self);
        rhs
    }
}

#[cfg(test)]
mod tests {
    use num_traits::Zero;
    use yui_polynomial::Poly2;

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
                TngComp::arc(2,4),
                TngComp::circ(1),
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

    #[test]
    fn stack_closed() {
        let mut c0 = Cob::from(CobComp::closed(0));
        let c1 = Cob::from(CobComp::closed(1));
        
        c0.stack(c1);

        assert_eq!(c0, Cob::new(vec![
            CobComp::closed(0),
            CobComp::closed(1)
        ]));
    }
    
    #[test]
    fn stack_cup_cap() {
        let mut c0 = Cob::from(CobComp::cup(TngComp::circ(0)));
        let c1 = Cob::from(CobComp::cap(TngComp::circ(0)));
        
        c0.stack(c1);

        assert_eq!(c0, Cob::new(vec![
            CobComp::sphere()
        ]));
    }
   
    #[test]
    fn stack_cap_cup() {
        let mut c0 = Cob::from(CobComp::cap(TngComp::circ(0)));
        let c1 = Cob::from(CobComp::cup(TngComp::circ(0)));
        
        c0.stack(c1);
        c0.sort_comps();

        assert_eq!(c0, Cob::new(vec![
            CobComp::cup(TngComp::circ(0)),
            CobComp::cap(TngComp::circ(0))
        ]));
    }
   
    #[test]
    fn stack_comps() {
        let mut c0 = Cob::new(vec![
            CobComp::id(TngComp::arc(0, 1)),
            CobComp::cup(TngComp::circ(2))
        ]);
        let c1 = Cob::new(vec![
            CobComp::cap(TngComp::circ(2)),
            CobComp::id(TngComp::arc(0, 1))
        ]);
        
        c0.stack(c1);
        c0.sort_comps();

        assert_eq!(c0, Cob::new(vec![
            CobComp::sphere(),
            CobComp::id(TngComp::arc(0, 1)),
        ]));
    }

    #[test]
    fn stack_id() {
        let c1 = Cob::new(vec![
            CobComp::from(&Crossing::from_pd_code([0,1,2,3])),
            CobComp::cup(TngComp::circ(4)),
            CobComp::cap(TngComp::circ(5)),
        ]);
        let c0 = Cob::id(&c1.src());
        let c2 = Cob::id(&c1.tgt());

        let mut e = c1.clone();
        e.stack(c2);

        assert_eq!(e, c1);

        let mut e = c0.clone();
        e.stack(c1.clone());
        assert_eq!(e, c1);
    }
   
    #[test]
    fn stack_torus() {
        let c0 = Cob::from(CobComp::cup(TngComp::Circ(0)));
        let c1 = Cob::new(vec![
            CobComp::split(
                TngComp::Circ(0),
                (TngComp::circ(1), TngComp::circ(2))
            )
        ]);
        let c2 = Cob::new(vec![
            CobComp::merge(
                (TngComp::circ(1), TngComp::circ(2)), 
                TngComp::Circ(3)
            )
        ]);
        let c3 = Cob::from(CobComp::cap(TngComp::Circ(3)));

        let mut c =  Cob::empty();
        c.stack(c0);

        assert_eq!(c.ncomps(), 1);
        assert_eq!(c.comp(0).src.ncomps(), 0);
        assert_eq!(c.comp(0).tgt.ncomps(), 1);
        assert_eq!(c.comp(0).genus, 0);

        c.stack(c1);
        assert_eq!(c.ncomps(), 1);
        assert_eq!(c.comp(0).src.ncomps(), 0);
        assert_eq!(c.comp(0).tgt.ncomps(), 2);
        assert_eq!(c.comp(0).genus, 0);

        c.stack(c2);
        assert_eq!(c.ncomps(), 1);
        assert_eq!(c.comp(0).src.ncomps(), 0);
        assert_eq!(c.comp(0).tgt.ncomps(), 1);
        assert_eq!(c.comp(0).genus, 1);

        c.stack(c3);
        assert_eq!(c.ncomps(), 1);
        assert_eq!(c.comp(0).src.ncomps(), 0);
        assert_eq!(c.comp(0).tgt.ncomps(), 0);
        assert_eq!(c.comp(0).genus, 1);
    }

    #[test]
    fn eval() { 
        type R = Poly2<'H', 'T', i32>;
        let h = R::variable(0);
        let t = R::variable(1);

        let c = CobComp::closed(0);
        assert_eq!(c.eval(&h, &t), R::zero());

        let c = CobComp::closed(1);
        assert_eq!(c.eval(&h, &t), R::from_const(2));

        let c = CobComp::closed(2);
        assert_eq!(c.eval(&h, &t), R::zero());

        let c = CobComp::closed(3);
        assert_eq!(c.eval(&h, &t), R::from_deg2_iter([((2, 0), 2), ((0, 1), 8)])); // 2(H^2 + 4T)
    }
}