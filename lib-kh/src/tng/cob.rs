//! Cobordism morphisms in Bar-Natan's category `Cob³_{/l}`: dotted surfaces
//! between Temperley–Lieb diagrams, modulo the local relations `S`, `T`, `4Tu`
//! plus the dotted skein `X² = h·X + t`. The second dot `Y = X − h` satisfies
//! `Y² = −h·Y + t` and `X·Y = t`. [`Cob`] is a connected-component
//! decomposition of a cobordism, and [`LcCob`] is its `R`-linear closure used
//! as the differential of [`super::TngComplex`].
//!
//! Reference:
//! - D. Bar-Natan, "Khovanov's homology for tangles and cobordisms",
//!   Geom. Topol. 9 (2005), 1443–1499.
//!   <https://doi.org/10.2140/gt.2005.9.1443>, <https://arxiv.org/abs/math/0410495>

use core::panic;
use std::fmt::Display;
use std::hash::Hash;
use std::collections::HashSet;
use std::ops::{Mul, MulAssign};
use auto_impl_ops::auto_ops;
use itertools::Itertools;
use num_traits::Zero;
use cartesian::cartesian;
use yui_core::util::format::subscript;
use yui_core::{AddMon, MathType, Ring, RingOps};
use yui_core::lc::{LcKey, Lc};
use yui_core::poly::Var2;
use super::tng::{Tng, TngComp};

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, derive_more::Display)]
pub enum Dot { 
    None, X, Y
}

#[derive(Clone, Copy, PartialEq, Eq, Debug, derive_more::Display)]
pub enum End { 
    Src, Tgt
}

#[derive(Clone, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub struct CobComp {
    src: Tng,
    tgt: Tng,
    genus: usize,
    dots: (usize, usize), // nums of X and Y dots resp.
    nb: usize,            // #∂-components — derived from (src, tgt)
}

impl CobComp { 
    fn new(src: Tng, tgt: Tng, genus: usize, dots: (usize, usize)) -> Self {
        let nb = Self::count_boundaries(&src, &tgt);
        Self::new_with_nb(src, tgt, genus, dots, nb)
    }

    fn new_with_nb(src: Tng, tgt: Tng, genus: usize, dots: (usize, usize), nb: usize) -> Self {
        debug_assert_eq!(src.end_pts().collect::<HashSet<_>>(), tgt.end_pts().collect());
        debug_assert_eq!(nb, Self::count_boundaries(&src, &tgt));
        Self { src, tgt, genus, dots, nb }
    }

    fn count_boundaries(src: &Tng, tgt: &Tng) -> usize {
        // Build a bitmask of arc-component indices in `t`, plus the count of
        // closed circles. Asserts `t.n_comps() <= 64`.
        let make = |t: &Tng| -> (u64, usize) {
            let n = t.n_comps();
            assert!(n <= 64, "count_boundaries: n_comps {n} exceeds u64 mask width");
            (0..n).fold((0u64, 0usize), |(arcs, circs), i| {
                if t.comp(i).is_arc() {
                    (arcs | (1u64 << i), circs)
                } else {
                    (arcs, circs + 1)
                }
            })
        };

        // First index `i` set in `mask` whose component in `t` is connectable to `c`.
        let next = |t: &Tng, mask: u64, c: &TngComp| -> Option<usize> {
            let mut m = mask;
            while m != 0 {
                let i = m.trailing_zeros() as usize;
                if t.comp(i).is_connectable(c) { return Some(i) }
                m &= m - 1;
            }
            None
        };

        let (mut src_arcs, src_circs) = make(src);
        let (mut tgt_arcs, tgt_circs) = make(tgt);

        debug_assert_eq!(src_arcs.count_ones(), tgt_arcs.count_ones());

        let mut side_circs = 0;
        while src_arcs != 0 {
            let mut i0 = src_arcs.trailing_zeros() as usize;
            loop {
                src_arcs &= !(1u64 << i0);
                let c0 = src.comp(i0);

                let j = next(tgt, tgt_arcs, c0).expect("no connectable tgt arc");
                tgt_arcs &= !(1u64 << j);
                let c1 = tgt.comp(j);

                match next(src, src_arcs, c1) {
                    Some(i1) => i0 = i1,
                    None => { side_circs += 1; break }
                }
            }
        }

        debug_assert_eq!(tgt_arcs, 0);

        src_circs + tgt_circs + side_circs
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

    pub fn cup(c: TngComp) -> Self {
        assert!(c.is_circle());
        Self::plain(
            Tng::empty(),
            Tng::from(c),
        )
    }

    pub fn cap(c: TngComp) -> Self {
        assert!(c.is_circle());
        Self::plain(
            Tng::from(c),
            Tng::empty(),
        )
    }

    pub fn with_dots(mut self, x: usize, y: usize) -> Self {
        self.dots = (x, y);
        self
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

    pub fn dots(&self) -> (usize, usize) { 
        self.dots
    }

    pub fn total_dots(&self) -> usize { 
        self.dots.0 + self.dots.1
    }

    pub fn n_boundaries(&self) -> usize {
        self.nb
    }

    pub fn end(&self, b: End) -> &Tng {
        match b {
            End::Src => &self.src,
            End::Tgt => &self.tgt
        }
    }

    pub fn end_mut(&mut self, b: End) -> &mut Tng {
        match b {
            End::Src => &mut self.src,
            End::Tgt => &mut self.tgt
        }
    }

    pub fn is_plain(&self) -> bool {
        self.dots == (0, 0) && self.genus == 0
    }

    pub fn is_closed(&self) -> bool {
        self.src.is_empty() &&
        self.tgt.is_empty()
    }

    pub fn is_zero_cob(&self) -> bool {
        self.is_closed() && 
        self.genus % 2 == 0 &&
        self.dots.0 == self.dots.1 // XY = T
    }

    pub fn is_removable(&self) -> bool {
        self.is_closed() &&
        self.genus == 0 &&
        (self.dots == (1, 0) || // ε.X.ι = 1,
         self.dots == (0, 1))   // ε.Y.ι = 1.
    }

    pub fn is_invertible(&self) -> bool {
        self.src.n_comps() == 1 &&
        self.tgt.n_comps() == 1 &&
        self.is_plain()
    }

    pub fn is_sdl(&self) -> bool {
        self.src.n_comps() == 2 && 
        self.tgt.n_comps() == 2 && 
        self.src.comps().all(|c| c.is_arc()) && 
        self.tgt.comps().all(|c| c.is_arc()) && 
        self.src != self.tgt && 
        self.genus == 0
    }

    pub fn inv(&self) -> Option<Self> { 
        if self.is_invertible() { 
            let inv = Self::plain(
                self.tgt.clone(),
                self.src.clone(),
            );
            Some(inv)
        } else {
            None
        }
    }

    // χ(S) = 2 - 2g(S) - #(∂S)
    pub fn euler_num(&self) -> i32 {
        let b = self.nb as i32;
        let g = self.genus as i32;
        2 - 2 * g - b
    }

    pub fn deg(&self) -> i32 { 
        let x = self.euler_num();
        let b = self.src.end_pts().count() as i32;
        let d = self.total_dots() as i32;
        x - (b / 2) - 2 * d
    }

    pub fn add_dot(&mut self, dot: Dot) { 
        match dot { 
            Dot::X => self.dots.0 += 1,
            Dot::Y => self.dots.1 += 1,
            _      => ()
        }
    }

    pub fn cap_off(&mut self, b: End, i: usize) {
        assert!(self.end(b).comp(i).is_circle());
        self.end_mut(b).remove_at(i);
        self.nb -= 1;
    }

    // connect = horizontal composition
    pub fn is_connectable(&self, other: &Self) -> bool { 
        self.src.comps().any(|c1| 
            if c1.is_arc() { 
                other.src.comps().any(|c2| { 
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

        let a = self.src.end_pts().filter(|e|
            other.src.end_pts().contains(&e)
        ).count() as i32;

        assert!(a > 0);

        let CobComp{ src, tgt, dots, .. } = other;

        self.src.connect(src);
        self.tgt.connect(tgt);
        self.nb = Self::count_boundaries(&self.src, &self.tgt);

        let b = self.nb as i32;
        let g = 2 - (x1 + x2 + b) + a;

        assert!(g >= 0);
        assert!(g % 2 == 0);

        self.genus = (g / 2) as usize;

        self.dots.0 += dots.0;
        self.dots.1 += dots.1;
    }

    pub fn should_reduce(&self) -> bool {
        self.is_zero_cob() ||
        self.is_removable() ||
        self.genus > 0 || 
        self.dots.0 >= 1 && self.dots.1 >= 1 ||
        self.dots.0 >= 2 ||
        self.dots.1 >= 2
    }

    pub fn reduce<R>(&self, h: &R, t: &R) -> LcCob<R>
    where R: Ring, for<'x> &'x R: RingOps<R> {
        fn eval<R>(c: &CobComp, nc: bool, g: usize, x: usize, y: usize, h: &R, t: &R) -> LcCob<R>
        where R: Ring, for<'x> &'x R: RingOps<R> {
            match (g, x, y) {
                // neck-cut (only valid when boundary has at most one component)
                (g, _, _) if g > 0 && nc => {
                    eval(c, nc, g-1, x+1, y, h, t) +
                    eval(c, nc, g-1, x, y+1, h, t)
                }

                // XY = t
                (0, x, y) if x >= 1 && y >= 1 =>
                    eval(c, nc, 0, x-1, y-1, h, t) * t,

                // X^2 = hX + t
                (0, x, 0) if x >= 2 =>
                    eval(c, nc, 0, x-1, 0, h, t) * h +
                    eval(c, nc, 0, x-2, 0, h, t) * t,

                // Y^2 = -hY + t
                (0, 0, y) if y >= 2 =>
                    eval(c, nc, 0, 0, y-1, h, t) * -h +
                    eval(c, nc, 0, 0, y-2, h, t) *  t,

                // XS = YS = 1
                (0, 1, 0) | (0, 0, 1) if c.is_closed() =>
                    Lc::from(Cob::empty()),

                // S = 0
                (0, 0, 0) if c.is_closed() =>
                    Lc::zero(),

                // default
                _ => {
                    let new = CobComp::new_with_nb(
                        c.src.clone(),
                        c.tgt.clone(),
                        g,
                        (x, y),
                        c.nb,
                    );
                    Lc::from(Cob::from(new))
                }
            }
        }

        let g = self.genus;
        let (x, y) = self.dots;
        let can_neck_cut = self.n_boundaries() <= 1;

        eval(self, can_neck_cut, g, x, y, h, t)
    }

    pub fn eval<R>(&self, h: &R, t: &R) -> R
    where R: Ring, for<'x> &'x R: RingOps<R> {
        assert!(self.is_closed(), "cannot eval: {}", self);

        let eval = self.reduce(h, t);

        assert!(eval.nterms() <= 1);

        if let Some((c, r)) = eval.any_term() { 
            assert!(c.is_empty());
            r.clone()
        } else { 
            R::zero()
        }
    }
}

impl Display for CobComp {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let dots = if self.total_dots() == 0 {
            String::new()
        } else {
            let (p, q) = self.dots;
            format!("{}・", Var2::<'X','Y', _>::from((p, q)))
        };


        if self.is_closed() { 
            return if self.genus == 0 { 
                write!(f, "{dots}S")
            } else { 
                write!(f, "{dots}Σ{}", subscript(self.genus as isize))                
            }
        }
        
        let base = match (self.src.n_comps(), self.tgt.n_comps()) {
            (0, 1) => "∪",
            (1, 0) => "∩",
            (1, 1) if self.is_invertible() => "I",
            (2, 1) => "∇",
            (1, 2) => "Δ",
            (2, 2) if self.is_sdl() => "X",
            _ => "Cob"
        }.to_string();
        
        let g = if self.genus == 0 { 
            "".to_string()
        } else { 
            format!(", g: {}", self.genus)
        };

        write!(f, "{dots}{base}({} -> {}{g})", self.src, self.tgt)
    }
}

#[derive(Clone, PartialEq, Eq, Hash, Debug, Default)]
pub struct Cob { 
    comps: Vec<CobComp>
}

impl Cob {
    pub fn new<I>(comps: I) -> Self
    where I: IntoIterator<Item = CobComp> { 
        let comps = comps.into_iter().sorted().collect_vec();
        Self { comps }
    }

    pub fn empty() -> Self { 
        Self::new(vec![])
    }
    
    pub fn id(v: &Tng) -> Self { 
        let comps = (0..v.n_comps()).map(|i| {
            let c = v.comp(i).clone();
            CobComp::id(c)
        });
        Self::new(comps)
    }

    pub fn n_comps(&self) -> usize { 
        self.comps.len()
    }

    pub fn comp(&self, i: usize) -> &CobComp { 
        &self.comps[i]
    }

    pub fn comp_mut(&mut self, i: usize) -> &mut CobComp { 
        &mut self.comps[i]
    }

    pub fn comps(&self) -> impl Iterator<Item = &CobComp> { 
        self.comps.iter()
    }

    pub fn find_comp(&mut self, b: End, c: &TngComp) -> Option<(usize, usize)> {
        self.comps.iter().enumerate().filter_map(|(i, comp)|
            comp.end(b).index_of(c).map(|p| (i, p))
        ).next()
    }

    pub fn n_boundaries(&self) -> usize { 
        self.comps.iter().map(|c| c.n_boundaries()).sum()
    }

    pub fn is_empty(&self) -> bool { 
        self.comps.is_empty()
    }

    pub fn is_zero_cob(&self) -> bool { 
        self.comps.iter().any(|c| c.is_zero_cob())
    }

    pub fn is_closed(&self) -> bool { 
        self.comps.iter().all(|c| c.is_closed())
    }

    pub fn is_invertible(&self) -> bool { 
        self.comps.iter().all(|c| c.is_invertible())
    }

    pub fn inv(&self) -> Option<Self> { 
        if self.is_invertible() { 
            let comps = self.comps.iter().map(|c| c.inv().unwrap());
            let inv = Self::new(comps);
            Some(inv)
        } else { 
            None
        }
    }

    pub fn euler_num(&self) -> i32 { 
        self.comps.iter().map(|c| c.euler_num()).sum()
    }

    pub fn deg(&self) -> i32 { 
        self.comps.iter().map(|c| c.deg()).sum()
    }

    pub fn cap_off(&mut self, b: End, c: &TngComp, x: Dot) {
        assert!(c.is_circle());
        let Some((i, p)) = self.find_comp(b, c) else { 
            panic!("{c} not found in {} ({b})", self)
        };

        let comp = self.comp_mut(i);
        comp.cap_off(b, p);
        comp.add_dot(x);

        if comp.is_removable() { 
            self.comps.remove(i);
        }

        self.normalize();
    }

    pub fn connect(&mut self, other: Cob) { // horizontal composition
        if other.is_empty() { return; }
        if self.is_empty() { *self = other; return; }

        let mut comps = std::mem::take(&mut self.comps);
        comps.reserve(other.comps.len());
        comps.extend(other.comps);

        while let Some(c) = Self::connect_next(&mut comps) {
            self.comps.push(c)
        }

        self.normalize()
    }

    // Take the next connected group from `comps` (transitively via shared
    // boundary endpoints) and return its merged CobComp. In-place frontier
    // index on `comps` — no auxiliary queue/Vec:
    //   comps[..unproc] = unprocessed
    //   comps[unproc..] = pending (in current group)
    // The pending region is drained back to empty before return.
    fn connect_next(comps: &mut Vec<CobComp>) -> Option<CobComp> {
        // Pop the seed and look for direct connectables. If none, the seed
        // is a singleton group and we can short-circuit (no genus rebuild).
        let seed = comps.pop()?;
        let mut unproc = comps.len();

        let mut i = 0;
        while i < unproc {
            if seed.is_connectable(&comps[i]) {
                comps.swap(i, unproc - 1);
                unproc -= 1;
            } else {
                i += 1;
            }
        }

        if comps.len() == unproc {
            return Some(seed);
        }

        // Real group: initialize accumulators from the seed.
        let mut dots = seed.dots;
        let mut x = seed.euler_num();
        let mut a = 0;
        let mut src = seed.src;
        let mut tgt = seed.tgt;

        while comps.len() > unproc {
            let cob = comps.pop().unwrap();

            let mut i = 0;
            while i < unproc {
                if cob.is_connectable(&comps[i]) {
                    comps.swap(i, unproc - 1);
                    unproc -= 1;
                } else {
                    i += 1;
                }
            }

            dots.0 += cob.dots.0;
            dots.1 += cob.dots.1;
            x += cob.euler_num();
            a += cob.src.end_pts().filter(|&e|
                src.end_pts().contains(&e)
            ).count();
            src.connect(cob.src);
            tgt.connect(cob.tgt);
        }

        let a = a as i32;
        let b = CobComp::count_boundaries(&src, &tgt) as i32;
        let g = 2 - (x + b) + a;

        assert!(g >= 0);
        assert!(g % 2 == 0);

        let genus = (g / 2) as usize;
        Some(CobComp::new_with_nb(src, tgt, genus, dots, b as usize))
    }

    pub fn is_stackable(&self, other: &Self) -> bool {
        self.comps.iter().fold(0, |n, c| n + c.tgt.n_comps()) ==
        other.comps.iter().fold(0, |n, c| n + c.src.n_comps()) &&
        self.comps.iter().all(|c| c.tgt.comps().all(|a|
            other.comps.iter().any(|c| c.src.contains(a))
        ))
    }

    pub fn stack(&mut self, other: Cob) { // vertical composition
        debug_assert!(
            self.is_stackable(&other),
            "{} cannot be stacked on {}", other, self
        );

        if self.is_empty() { 
            *self = other;
            return;
        } else if other.is_empty() { 
            return;
        }

        let mut bot = std::mem::take(&mut self.comps);
        let mut top = other.comps;

        while let Some(c) = Self::stack_next(&mut bot, &mut top) {
            self.comps.push(c)
        }

        self.normalize()
    }

    // Take the next connected group from `bot` & `top` and return its merged
    // CobComp. Uses an in-place frontier index instead of an explicit queue:
    //   bot[..bot_unproc] = unprocessed (available for later groups)
    //   bot[bot_unproc..] = pending (in current group, awaiting absorption)
    // Same for top. Both pending regions are empty again on return.
    fn stack_next(bot: &mut Vec<CobComp>, top: &mut Vec<CobComp>) -> Option<CobComp> {
        if bot.is_empty() && top.is_empty() { return None }

        let mut src = Tng::empty();
        let mut tgt = Tng::empty();
        let mut dots = (0, 0);
        let mut x = 0 as i32;
        let mut a = 0 as i32;

        let mut bot_unproc = bot.len();
        let mut top_unproc = top.len();

        // Seed: trailing comp becomes the first pending member of the group.
        if bot_unproc > 0 {
            bot_unproc -= 1;
        } else {
            top_unproc -= 1;
        }

        while bot.len() > bot_unproc || top.len() > top_unproc {
            if bot.len() > bot_unproc {
                let cob = bot.pop().unwrap();
                for c in cob.tgt.comps() {
                    if let Some(i) = top[..top_unproc].iter().position(|t| t.src.contains(c)) {
                        top.swap(i, top_unproc - 1);
                        top_unproc -= 1;
                    }
                    if c.is_arc() {
                        a += 1;
                    }
                }
                dots.0 += cob.dots.0;
                dots.1 += cob.dots.1;
                x += cob.euler_num();
                src.connect(cob.src);
            } else {
                let cob = top.pop().unwrap();
                for c in cob.src.comps() {
                    if let Some(i) = bot[..bot_unproc].iter().position(|b| b.tgt.contains(c)) {
                        bot.swap(i, bot_unproc - 1);
                        bot_unproc -= 1;
                    }
                }
                dots.0 += cob.dots.0;
                dots.1 += cob.dots.1;
                x += cob.euler_num();
                tgt.connect(cob.tgt);
            }
        }

        let b = CobComp::count_boundaries(&src, &tgt) as i32;
        let g = 2 - (x + b) + a;

        assert!(g >= 0);
        assert!(g % 2 == 0);

        let genus = (g / 2) as usize;
        Some(CobComp::new_with_nb(src, tgt, genus, dots, b as usize))
    }

    pub fn should_reduce(&self) -> bool {
        self.comps.iter().any(|c| c.should_reduce())
    }

    pub fn reduce<R>(mut self, h: &R, t: &R) -> Lc<Cob, R>
    where R: Ring, for<'x> &'x R: RingOps<R> {
        if self.is_zero_cob() {
            return Lc::zero()
        }
        if !self.should_reduce() {
            return Lc::from(self)
        }
        if self.comps.len() == 1 {
            return self.comps.into_iter().next().unwrap().reduce(h, t);
        }

        let need_reduce: Vec<_> = self.comps.extract_if(.., |c| c.should_reduce()).collect();
        let init = LcCob::from(self);
        
        need_reduce.into_iter().fold(init, |res, c| {
            let e = c.reduce(h, t);
            debug_assert!(e.keys().all(|c| c.n_comps() <= 1));

            res.apply_bilin(&e, |c1, c2| {
                let mut c = c1.clone();
                if let Some(c2) = c2.comps.first() {
                    c.comps.push(c2.clone());
                }
                c
            })
        }).map_keys(|mut c| {
            c.normalize();
            c
        })
    }

    pub fn eval<R>(&self, h: &R, t: &R) -> R
    where R: Ring, for<'x> &'x R: RingOps<R> {
        let comps = self.comps.iter().map(|c| 
            c.eval(h, t)
        );
        R::product(comps)
    }

    fn normalize(&mut self) {
        self.comps.retain(|c| !c.is_removable());
        self.comps.sort()
    }

    #[cfg(debug_assertions)]
    pub fn reconst_src(&self) -> Tng { 
        self.comps.iter().fold(Tng::empty(), |mut t, c| {
            t.connect(c.src.clone());
            t
        })
    }

    #[cfg(debug_assertions)]
    pub fn reconst_tgt(&self) -> Tng { 
        self.comps.iter().fold(Tng::empty(), |mut t, c| {
            t.connect(c.tgt.clone());
            t
        })
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
            write!(f, "(∅)")
        } else if self.comps.len() == 1 {
            write!(f, "{}", self.comps[0])
        } else { 
            let cobs = self.comps.iter().join(" ⊔ ");
            write!(f, "|{cobs}|")
        }
    }
}

impl PartialOrd for Cob {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Cob {
    fn cmp(&self, _other: &Self) -> std::cmp::Ordering {
        // TODO
        std::cmp::Ordering::Equal
    }
}

impl MathType for Cob {
    fn math_symbol() -> String {
        "Cob".to_string()
    }
}

impl LcKey for Cob {}

#[auto_ops]
impl Mul for Cob {
    type Output = Cob;
    fn mul(self, mut rhs: Self) -> Self::Output {
        rhs.stack(self);
        rhs
    }
}

pub type LcCob<R> = Lc<Cob, R>; // R-linear combination of cobordisms.

pub trait LcCobTrait: Sized {
    type R;
    fn is_closed(&self) -> bool;
    fn is_invertible(&self) -> bool;
    fn is_stackable(&self, other: &Self) -> bool;
    fn inv(&self) -> Option<Self>;
    fn connect(self, c: &Cob) -> Self;
    fn connect_ref(&self, c: &Cob) -> Self;
    fn cap_off(self, b: End, c: &TngComp, dot: Dot) -> Self;
    fn should_reduce(&self) -> bool;
    fn reduce(self, h: &Self::R, t: &Self::R) -> Self;
    fn eval(&self, h: &Self::R, t: &Self::R) -> Self::R;
}

impl<R> LcCobTrait for LcCob<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;

    fn is_closed(&self) -> bool { 
        self.iter().all(|(f, _)| f.is_closed())
    }

    fn is_invertible(&self) -> bool { 
        self.nterms() == 1 && 
        self.iter().next().map(|(c, a)| 
            c.is_invertible() && a.is_unit()
        ).unwrap_or(false)
    }

    fn is_stackable(&self, other: &Self) -> bool { 
        cartesian!(self.keys(), other.keys()).all(|(a, b)| 
            a.is_stackable(b)
        )
    }

    fn inv(&self) -> Option<Self> { 
        if let Some((Some(cinv), Some(ainv))) = self.iter().next().map(|(c, a)| 
            (c.inv(), a.inv())
        ) { 
            let inv = LcCob::from((cinv, ainv));
            Some(inv)
        } else { 
            None
        }
    }

    fn connect(self, c: &Cob) -> Self {
        mut_cob(self, |cob| cob.connect(c.clone()))
    }

    fn connect_ref(&self, c: &Cob) -> Self {
        mut_cob_ref(self, |cob| cob.connect(c.clone()))
    }

    fn cap_off(self, b: End, c: &TngComp, dot: Dot) -> Self {
        mut_cob(self, |cob| cob.cap_off(b, c, dot) )
    }

    fn should_reduce(&self) -> bool {
        self.keys().any(|c| c.should_reduce())
    }

    fn reduce(self, h: &Self::R, t: &Self::R) -> Self {
        if self.should_reduce() { 
            LcCob::sum(self.into_iter().map(|(cob, r)|
                cob.reduce(h, t) * r
            ))
        } else { 
            self
        }
    }

    fn eval(&self, h: &R, t: &R) -> R {
        let coeffs = self.iter().map(|(c, a)|
            a * c.eval(h, t)
        );
        R::sum(coeffs)
    }
}

/// Apply `f` to each `Cob` in `this`, dropping any terms that become zero.
fn mut_cob<R, F>(this: LcCob<R>, f: F) -> LcCob<R>
where
    R: Ring, for<'x> &'x R: RingOps<R>,
    F: Fn(&mut Cob),
{
    this.into_iter().filter_map(|(mut cob, r)| {
        f(&mut cob);
        (!cob.is_zero_cob()).then_some((cob, r))
    }).collect()
}

/// Non-consuming variant of [`modify_cob`]: clones each `Cob` individually
/// rather than the whole `LcCob`.
fn mut_cob_ref<R, F>(this: &LcCob<R>, f: F) -> LcCob<R>
where
    R: Ring, for<'x> &'x R: RingOps<R>,
    F: Fn(&mut Cob),
{
    this.iter().filter_map(|(cob, r)| {
        let mut new_cob = cob.clone();
        f(&mut new_cob);
        (!new_cob.is_zero_cob()).then_some((new_cob, r.clone()))
    }).collect()
}

#[cfg(test)]
mod tests {
    use num_traits::Zero;
    use maplit::hashmap;
    use yui_core::CloneAnd;
    use yui_core::poly::Poly2;
    use yui_core::bitseq::Bit;
    use yui_link::Node;

    use super::CobComp;
    use super::*;

    fn sdl(r0: (TngComp, TngComp), r1: (TngComp, TngComp)) -> CobComp {
        CobComp::plain(
            Tng::new(vec![r0.0, r0.1]),
            Tng::new(vec![r1.0, r1.1]),
        )
    }

    fn pants(from: (TngComp, TngComp), to: TngComp) -> CobComp {
        assert!(from.0.is_circle() || from.1.is_circle());
        CobComp::plain(
            Tng::new(vec![from.0, from.1]),
            Tng::new(vec![to]),
        )
    }

    fn copants(from: TngComp, to: (TngComp, TngComp)) -> CobComp {
        assert!(to.0.is_circle() || to.1.is_circle());
        CobComp::plain(
            Tng::new(vec![from]),
            Tng::new(vec![to.0, to.1]),
        )
    }

    fn closed(g: usize) -> CobComp {
        CobComp::new(Tng::empty(), Tng::empty(), g, (0, 0))
    }

    #[test]
    fn is_connectable() {
        let src = Tng::new(vec![
            TngComp::arc([1, 2]),
            TngComp::arc([3, 4]),
            TngComp::circ([10]),
        ]);
        let tgt = Tng::new(vec![
            TngComp::arc([1, 3]),
            TngComp::arc([2, 4]),
            TngComp::circ([11]),
        ]);
        let c = CobComp::plain(src, tgt);

        let c1 = CobComp::id(
            TngComp::arc([0, 1])
        );
        let c2 = sdl(
            (TngComp::arc([0, 1]), TngComp::arc([90, 91])),
            (TngComp::arc([0, 90]), TngComp::arc([1, 91])),
        );
        let c3 = CobComp::id(TngComp::arc([5, 6]));

        assert!(c.is_connectable(&c1));
        assert!(c.is_connectable(&c2));
        assert!(!c.is_connectable(&c3));
    }

    #[test]
    fn connect1() { 
        let src = Tng::new(vec![
            TngComp::arc([1, 2]),
            TngComp::arc([3, 4]),
            TngComp::circ([10]),
        ]);
        let tgt = Tng::new(vec![
            TngComp::arc([1, 3]),
            TngComp::arc([2, 4]),
            TngComp::circ([11]),
        ]);

        let mut c = CobComp::plain(src, tgt);
        c.connect(CobComp::id(
            TngComp::arc([0, 1])
        ));

        assert_eq!(c, CobComp::plain(
            Tng::new(vec![
                TngComp::arc([0, 1, 2]),
                TngComp::arc([3, 4]),
                TngComp::circ([10]),
            ]),
            Tng::new(vec![
                TngComp::arc([0, 1, 3]),
                TngComp::arc([2, 4]),
                TngComp::circ([11]),
            ]),
        ));
    }

    #[test]
    fn connect2() { 
        let src = Tng::new(vec![
            TngComp::arc([1, 2]),
            TngComp::arc([3, 4]),
            TngComp::circ([10]),
        ]);
        let tgt = Tng::new(vec![
            TngComp::arc([1, 3]),
            TngComp::arc([2, 4]),
            TngComp::circ([11]),
        ]);

        let mut c = CobComp::plain(src, tgt);
        c.connect(CobComp::id(
            TngComp::arc([1, 3])
        ));

        assert_eq!(c, CobComp::plain(
            Tng::new(vec![
                TngComp::arc([2, 1, 3, 4]),
                TngComp::circ([10]),
            ]),
            Tng::new(vec![
                TngComp::arc([2, 4]),
                TngComp::circ([1, 3]),
                TngComp::circ([11]),
            ]),
        ));
    }

    #[test]
    fn euler_num() { 
        let c0 = CobComp::id(
            TngComp::arc([1, 2])
        );
        let c1 = sdl(
            (TngComp::arc([3, 4]), TngComp::arc([5, 6])),
            (TngComp::arc([4, 5]), TngComp::arc([6, 3])),
        );
        let c2 = CobComp::plain(
            Tng::from(TngComp::circ([10])),
            Tng::new(vec![TngComp::circ([10]), TngComp::circ([11])]),
        );
        let c3 = CobComp::cup(
            TngComp::circ([20])
        );
        let c4 = CobComp::cap(
            TngComp::circ([30])
        );

        assert_eq!(c0.n_boundaries(), 1);
        assert_eq!(c1.n_boundaries(), 1);
        assert_eq!(c2.n_boundaries(), 3);
        assert_eq!(c3.n_boundaries(), 1);
        assert_eq!(c4.n_boundaries(), 1);

        assert_eq!(c0.euler_num(), 1);
        assert_eq!(c1.euler_num(), 1);
        assert_eq!(c2.euler_num(), -1);
        assert_eq!(c3.euler_num(), 1);
        assert_eq!(c4.euler_num(), 1);

        let cob = Cob::new(vec![c0,c1,c2,c3,c4]);
        assert_eq!(cob.euler_num(), 3);
        assert_eq!(cob.n_boundaries(), 7);
    }

    #[test]
    fn connect_incr_genus() { 
        let mut c0 = CobComp::plain(
            Tng::new(vec![
                TngComp::arc([1, 2]),
                TngComp::arc([3, 4])
            ]),
            Tng::new(vec![
                TngComp::arc([1, 2]),
                TngComp::arc([3, 4])
            ]),
        );
        let c1 = CobComp::id(
            TngComp::arc([1, 3])
        );
        let c2 = CobComp::id(
            TngComp::arc([2, 4])
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

        c0.cap_off(End::Src, 0);

        assert_eq!(c0.genus, 1);
        assert_eq!(c0.euler_num(), -1);

        c0.cap_off(End::Tgt, 0);

        assert_eq!(c0.genus, 1);
        assert_eq!(c0.euler_num(), 0);
        assert!(c0.is_closed()); // torus
    }

    #[test]
    fn inv() { 
        let cc0 = CobComp::id(TngComp::arc([0, 1]));
        let cc1 = CobComp::plain(
            Tng::from(TngComp::circ([2])),
            Tng::from(TngComp::circ([3])),
        );

        assert!(cc0.is_invertible());
        assert_eq!(cc0.inv(), Some(cc0.clone()));

        assert!(cc1.is_invertible());
        assert_eq!(cc1.inv(), Some(CobComp::plain(
            Tng::from(TngComp::circ([3])),
            Tng::from(TngComp::circ([2])),
        )));

        let c0 = Cob::new(vec![cc0, cc1]);

        assert!(c0.is_invertible());
        assert_eq!(c0.inv(), Some(Cob::new(vec![
            c0.comp(0).inv().unwrap(),
            c0.comp(1).inv().unwrap()
        ])));

        let c1 = Cob::from(
            sdl(
                (TngComp::arc([1, 2]), TngComp::arc([3, 4])),
                (TngComp::arc([1, 3]), TngComp::arc([2, 4])),
            )
        );

        assert!(!c1.is_invertible());
        assert_eq!(c1.inv(), None);

        let c2 = c0.clone_and(|c2|
            c2.comps[0].add_dot(Dot::X)
        );

        assert!(!c2.is_invertible());
        assert_eq!(c2.inv(), None);

        let c3 = c0.clone_and(|c3|
            c3.comps[0].genus += 1
        );

        assert!(!c3.is_invertible());
        assert_eq!(c3.inv(), None);
    }

    #[test]
    fn mor_inv() { 
        let c = Cob::id(&Tng::new(vec![
            TngComp::arc([0, 1]),
            TngComp::arc([2, 3])
        ]));
        let f = LcCob::from((c.clone(), -1));

        assert!(f.is_invertible());
        assert_eq!(f.inv(), Some(f.clone()));

        let f = LcCob::from((c.clone(), 2));
        assert!(!f.is_invertible());
        assert_eq!(f.inv(), None);
    }

    #[test]
    fn stack_closed() {
        let mut c0 = Cob::from(closed(0));
        let c1 = Cob::from(closed(1));
        
        c0.stack(c1);

        assert_eq!(c0, Cob::new(vec![
            closed(0),
            closed(1)
        ]));
    }
    
    #[test]
    fn stack_cup_cap() {
        let mut c0 = Cob::from(CobComp::cup(TngComp::circ([0])));
        let c1 = Cob::from(CobComp::cap(TngComp::circ([0])));
        
        c0.stack(c1);

        assert_eq!(c0, Cob::new(vec![
            closed(0)
        ]));
    }
   
    #[test]
    fn stack_cap_cup() {
        let mut c0 = Cob::from(CobComp::cap(TngComp::circ([0])));
        let c1 = Cob::from(CobComp::cup(TngComp::circ([0])));
        
        c0.stack(c1);

        assert_eq!(c0, Cob::new(vec![
            CobComp::cup(TngComp::circ([0])),
            CobComp::cap(TngComp::circ([0]))
        ]));
    }
   
    #[test]
    fn stack_comps() {
        let mut c0 = Cob::new(vec![
            CobComp::id(TngComp::arc([0, 1])),
            CobComp::cup(TngComp::circ([2]))
        ]);
        let c1 = Cob::new(vec![
            CobComp::cap(TngComp::circ([2])),
            CobComp::id(TngComp::arc([0, 1]))
        ]);
        
        c0.stack(c1);

        assert_eq!(c0, Cob::new(vec![
            closed(0),
            CobComp::id(TngComp::arc([0, 1])),
        ]));
    }

    #[test]
    fn stack_id() {
        let node = Node::from_pd_code([1,4,2,5]);
        let c1 = Cob::new(vec![
            CobComp::plain(
                Tng::from_resolved(&node.resolve(Bit::Bit0), None),
                Tng::from_resolved(&node.resolve(Bit::Bit1), None),
            ),
            CobComp::cup(TngComp::circ([10])),
            CobComp::cap(TngComp::circ([11])),
        ]);
        let c0 = Cob::id(&c1.reconst_src());
        let c2 = Cob::id(&c1.reconst_tgt());

        let e = c1.clone_and(|e|
            e.stack(c2)
        );

        assert_eq!(e, c1);

        let e = c0.clone_and(|e|
            e.stack(c1.clone())
        );
        assert_eq!(e, c1);
    }
   
    #[test]
    fn stack_torus() {
        let c0 = Cob::from(CobComp::cup(TngComp::circ([0])));
        let c1 = Cob::new(vec![
            copants(
                TngComp::circ([0]),
                (TngComp::circ([1]), TngComp::circ([2]))
            )
        ]);
        let c2 = Cob::new(vec![
            pants(
                (TngComp::circ([1]), TngComp::circ([2])),
                TngComp::circ([3])
            )
        ]);
        let c3 = Cob::from(CobComp::cap(TngComp::circ([3])));

        let mut c =  Cob::empty();
        c.stack(c0);

        assert_eq!(c.n_comps(), 1);
        assert_eq!(c.comp(0).src.n_comps(), 0);
        assert_eq!(c.comp(0).tgt.n_comps(), 1);
        assert_eq!(c.comp(0).genus, 0);

        c.stack(c1);
        assert_eq!(c.n_comps(), 1);
        assert_eq!(c.comp(0).src.n_comps(), 0);
        assert_eq!(c.comp(0).tgt.n_comps(), 2);
        assert_eq!(c.comp(0).genus, 0);

        c.stack(c2);
        assert_eq!(c.n_comps(), 1);
        assert_eq!(c.comp(0).src.n_comps(), 0);
        assert_eq!(c.comp(0).tgt.n_comps(), 1);
        assert_eq!(c.comp(0).genus, 1);

        c.stack(c3);
        assert_eq!(c.n_comps(), 1);
        assert_eq!(c.comp(0).src.n_comps(), 0);
        assert_eq!(c.comp(0).tgt.n_comps(), 0);
        assert_eq!(c.comp(0).genus, 1);
    }

    #[test]
    fn eval() { 
        type R = Poly2<'H', 'T', i32>;
        
        let ht = R::mono;
        let h = R::variable(0);
        let t = R::variable(1);

        let c = closed(0);
        assert_eq!(c.eval(&h, &t), R::zero());

        let c = closed(1);
        assert_eq!(c.eval(&h, &t), R::from_const(2));

        let c = closed(2);
        assert_eq!(c.eval(&h, &t), R::zero());

        let c = closed(3);
        assert_eq!(c.eval(&h, &t), R::from_iter([(ht(2, 0), 2), (ht(0, 1), 8)])); // 2(H^2 + 4T)
    }

    #[test]
    fn reduce() { 
        let mut c0 = CobComp::id(TngComp::circ([1]));
        c0.add_dot(Dot::X);

        let mut c1 = CobComp::id(TngComp::circ([1]));
        c1.add_dot(Dot::X);
        c1.add_dot(Dot::X);

        let c = LcCob::from_iter(hashmap! { 
            Cob::from(c0) => -2,
            Cob::from(c1) => 1
        });

        assert!(!c.is_zero());
        assert!(c.reduce(&2, &0).is_zero()); // X^2 = 2X
    }

    #[test]
    fn cobcomp_display() {
        // Closed shapes: "S" for genus-0, "Σ_n" for higher genus.
        assert_eq!(closed(0).to_string(),                "S");
        assert_eq!(closed(0).with_dots(1, 0).to_string(), "X・S");
        assert_eq!(closed(0).with_dots(0, 1).to_string(), "Y・S");
        assert_eq!(closed(0).with_dots(1, 1).to_string(), "XY・S");
        assert_eq!(closed(0).with_dots(2, 0).to_string(), "X²・S");
        assert_eq!(closed(1).to_string(),                "Σ₁");
        assert_eq!(closed(2).with_dots(1, 0).to_string(), "X・Σ₂");

        // Open shapes: base symbol + (src -> tgt) + optional ", g: N".
        let c = TngComp::circ([0]);
        assert_eq!(CobComp::cup(c.clone()).to_string(), "∪(∅ -> ⚪︎(0))");
        assert_eq!(CobComp::cap(c.clone()).to_string(), "∩(⚪︎(0) -> ∅)");
        assert_eq!(CobComp::id(c.clone()).to_string(),  "I(⚪︎(0) -> ⚪︎(0))");

        // Saddle.
        let sdl = CobComp::plain(
            Tng::new(vec![TngComp::arc([0, 1]), TngComp::arc([2, 3])]),
            Tng::new(vec![TngComp::arc([0, 2]), TngComp::arc([1, 3])]),
        );
        assert_eq!(sdl.to_string(), "X({[0-1], [2-3]} -> {[0-2], [1-3]})");

        // Genus on open cob: drops to the default `Cob` label; trailing ", g: N".
        let mut handle = CobComp::id(c);
        handle.genus = 1;
        assert_eq!(handle.to_string(), "Cob(⚪︎(0) -> ⚪︎(0), g: 1)");
    }
}