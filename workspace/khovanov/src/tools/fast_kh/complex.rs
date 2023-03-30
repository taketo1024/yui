// #![allow(unused)] // TODO remove 

use std::fmt::Display;
use std::ops::RangeInclusive;

use num_traits::Zero;
use cartesian::cartesian;
use yui_core::{Ring, RingOps};
use yui_homology::GenericChainComplex;
use yui_lin_comb::LinComb;
use yui_matrix::sparse::{MatType, SpMat};
use yui_matrix::dense::Mat;
use yui_link::{Crossing, Resolution, State, Component};

use crate::KhAlgLabel;
use crate::tools::fast_kh::cob::CobComp;

use super::cob::{Cob, Dot, End};
use super::tng::{Tng, TngUpdate};

#[derive(Clone)]
pub struct Obj { 
    state: State, 
    label: Vec<KhAlgLabel>,
    tangle: Tng,
}

impl Obj { 
    fn new(state: State, label: Vec<KhAlgLabel>, tangle: Tng) -> Self { 
        Self { state, label, tangle }
    }

    fn empty() -> Self { 
        Self::new(State::empty(), vec![], Tng::empty())
    }

    fn append_bit(&mut self, r: Resolution) { 
        self.state.append_b(r)
    }

    fn append_label(&mut self, x: KhAlgLabel) { 
        self.label.push(x)
    }

    fn append_arc(&mut self, a: Component) -> TngUpdate {
        self.tangle.append_arc(a)
    }

    fn find_loop(&self) -> Option<usize> { 
        self.tangle.find_loop()
    }

    fn deloop(&mut self, i: usize) {
        self.tangle.deloop(i);
    }
}

impl Display for Obj {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{{{}; {:?}; {}}}", self.state, self.label, self.tangle)
    }
}

type Mor = LinComb<Cob, i32>; // Z-linear combination of cobordisms.

trait MorTrait: Sized {
    fn is_invertible(&self) -> bool;
    fn inv(&self) -> Option<Self>;
    fn cap_off(self, r: usize, dot: Dot, e: End) -> Self;
    fn eval<R>(&self, h: &R, t: &R) -> R
    where R: Ring + From<i32>, for<'x> &'x R: RingOps<R>;
}

impl MorTrait for Mor {
    fn is_invertible(&self) -> bool { 
        self.len() == 1 && 
        self.iter().next().map(|(c, a)| 
            c.is_invertible() && a.is_unit()
        ).unwrap_or(false)
    }

    fn inv(&self) -> Option<Self> { 
        if let Some((Some(cinv), Some(ainv))) = self.iter().next().map(|(c, a)| 
            (c.inv(), a.inv())
        ) { 
            let inv = Mor::from((cinv, ainv));
            Some(inv)
        } else { 
            None
        }
    }

    fn cap_off(self, p: usize, dot: Dot, e: End) -> Self {
        self.into_map(|mut cob, r| { 
            cob.cap_off(p, dot, e);

            if cob.is_zero() { 
                (cob, 0)
            } else { 
                (cob, r)
            }
        })
    }

    fn eval<R>(&self, h: &R, t: &R) -> R
    where R: Ring + From<i32>, for<'x> &'x R: RingOps<R> {
        self.iter().map(|(c, &a)| { 
            R::from(a) * c.eval(h, t)
        }).sum()
    }
}
pub struct TngComplex {
    objs: Vec<Vec<Obj>>, // ⊕_i ⊕_s T[s]
    mats: Vec<Mat<Mor>>
}

impl TngComplex {
    pub fn new() -> Self { 
        let objs = vec![vec![Obj::empty()]];
        let mats = vec![Mat::zero((0, 1))];
        TngComplex{ objs, mats }
    }

    pub fn len(&self) -> usize { 
        self.objs.len()
    }

    pub fn rank(&self, i: usize) -> usize { 
        if i < self.objs.len() { 
            self.objs[i].len()
        } else { 
            0
        }
    }

    pub fn obj(&self, i: usize, j: usize) -> &Obj { 
        &self.objs[i][j]
    }

    pub fn tng(&self, i: usize, j: usize) -> &Tng { 
        &self.objs[i][j].tangle
    }

    pub fn mat(&self, i: usize) -> &Mat<Mor> { 
        &self.mats[i]
    }

    pub fn append(&mut self, x: &Crossing) {
        self.modify(|b| b.append(x));
    }

    pub fn deloop(&mut self) { 
        self.modify(|b| b.deloop());
    }

    fn modify<F>(&mut self, f: F)
    where F: Fn(&mut TngComplexBuilder) { 
        let objs = std::mem::take(&mut self.objs);
        let mats = std::mem::take(&mut self.mats);

        let mut builder = TngComplexBuilder::new(objs, mats);
        f(&mut builder);

        (self.objs, self.mats) = builder.result();
    }

    pub fn is_completely_delooped(&self) -> bool { 
        self.objs.iter().all(|c|
            c.iter().all(|v|
                v.tangle.is_empty()
            )
        )
    }

    pub fn as_generic<'a, R>(&'a self, h: &R, t: &R) -> GenericChainComplex<R, RangeInclusive<isize>> 
    where R: Ring + From<i32>, for<'x> &'x R: RingOps<R> {
        assert!(self.is_completely_delooped());

        let n = self.len();
        GenericChainComplex::ascending((0..n).map( |i| {
            let d = self.mat(i);
            SpMat::generate(d.shape(), |set| { 
                for (i, j, c) in d.iter() { 
                    if c.is_zero() { continue }
                    let a = c.eval(h, t);
                    set(i, j, a)
                }
            })
        }).collect())
    }

    pub fn describe(&self) { 
        let mut str = "".to_string();
        for (i, ci) in self.objs.iter().enumerate() { 
            str += &format!("C[{i}]:\n");
            for (j, v) in ci.iter().enumerate() { 
                str += &format!(" - [{j}]: {v}\n");
            }
            str += "\n";
        }
        str += "\n";
        for (i, di) in self.mats.iter().enumerate() { 
            str += &format!("d[{i}]:\n{}\n\n", di);
        }
        println!("{str}");
    }
}

struct TngComplexBuilder { 
    objs: Vec<Vec<Obj>>, // ⊕_i ⊕_s T[s]
    mats: Vec<Mat<Mor>>,
    ranks: Vec<usize>,
    updates: Vec<Vec<Vec<TngUpdate>>> // tmp data used during `append`.
}

impl TngComplexBuilder { 
    fn new(objs: Vec<Vec<Obj>>, mats: Vec<Mat<Mor>>) -> Self { 
        let ranks = vec![];
        let updates = vec![];
        Self { objs, mats, ranks, updates }
    }

    fn result(self) -> (Vec<Vec<Obj>>, Vec<Mat<Mor>>) { 
        (self.objs, self.mats)
    }

    //                          d0
    //                 C[i]#x0 ---> C[i+1]#x0
    //                   |             :
    //                 f |             :
    //            -d1    V             :
    //  C[i-1]#x1 ---> C[i]#x1 .... C[i+1]#x1 
    //
    //  C'[i] = C[i]#x0 ⊕ C[i-1]#x1
    //     d' = [d0    ]
    //          [f  -d1]

    fn append(&mut self, x: &Crossing) {
        assert!(self.updates.is_empty());

        self.init_updates();
        self.extend_objs();
        self.extend_mats();

        self.modif_objs(x);
        self.modif_mats();
    }

    fn init_updates(&mut self) { 
        self.ranks.clear();
        self.updates.clear();

        let n = self.objs.len();

        for i in 0 .. n { 
            let r = self.objs[i].len();
            self.ranks.push(r);
        }
        self.ranks.push(0);

        for i in 0 .. n+1 { 
            let r = if i > 0 { 
                self.ranks[i] + self.ranks[i - 1]
            } else { 
                self.ranks[0]
            };
            let slot = vec![vec![]; r];
            self.updates.push(slot);
        }
    }

    fn extend_objs(&mut self) {
        let n = self.objs.len();
        let mut clone = self.objs.clone();

        // duplicate objs with deg-shift.
        self.objs.push(vec![]);

        for i in 0..n { 
            self.objs[i+1].append(&mut clone[i]);
        }
    }

    fn extend_mats(&mut self) { 
        let n = self.mats.len();
        let mut mats = vec![];

        for i in 0..n { 
            let r = self.ranks[i];

            let d0 = &self.mats[i];
            let d1 = if i > 0 {  
                -&self.mats[i-1]
            } else { 
                Mat::zero((r, 0))
            };

            let f = Mat::diag((r, r), (0..r).map(|j| {
                let v = &self.objs[i][j];
                let f = Cob::id_for(&v.tangle);
                Mor::from(f)
            }));
            let o = Mat::zero((d0.rows(), d1.cols()));
            let d = Mat::combine_blocks([
                &d0, &o, 
                &f,  &d1
            ]);

            mats.push(d)
        }

        let d_n = Mat::zero((0, self.mats[n-1].cols()));
        mats.push(d_n);

        self.mats = mats;
    }

    fn modif_objs(&mut self, x: &Crossing) { 
        let n = self.objs.len() - 1;
        for i in 0..n { 
            let r = self.ranks[i];
            let s = self.ranks[i+1];
            for j in 0..r { 
                self.append_x_res(i,   j,   x, Resolution::Res0); // C[i]#x0
                self.append_x_res(i+1, s+j, x, Resolution::Res1); // C[i]#x1
            }
        }
    }

    fn append_x_res(&mut self, i: usize, j: usize, x: &Crossing, r: Resolution) {
        let v = &mut self.objs[i][j];
        v.append_bit(r);

        let (r0, r1) = x.res_arcs(r);
        self.append_arc(i, j, r0);
        self.append_arc(i, j, r1);
    }

    fn append_arc(&mut self, i: usize, j: usize, arc: Component) {
        let v = &mut self.objs[i][j];
        let e = v.append_arc(arc);
        self.updates[i][j].push(e);
    }

    fn modif_mats(&mut self) { 
        let n = self.mats.len() - 1;
        for i in 0..n { 
            let r = self.ranks[i];
            let s = self.ranks[i+1];

            // modify d0
            for (j, k) in cartesian!(0..r, 0..s) { 
                self.insert_cyl(i, j, k);
            }

            // modify f
            for j in 0..r { 
                self.insert_sdl(i, j, s+j);
            }

            // modify d1
            if i > 0 { 
                let q = self.ranks[i-1];
                for (j, k) in cartesian!(0..q, 0..r) { 
                    self.insert_cyl(i, r+j, s+k);
                }
            }
        }
    }

    fn insert_cyl(&mut self, i: usize, j: usize, k: usize) { 
        self.insert_cob(i, j, k, |r0, r1| { 
            let c0 = CobComp::cyl(r0.0, r1.0);
            let c1 = CobComp::cyl(r0.1, r1.1);
            vec![c0, c1]
        })
    }

    fn insert_sdl(&mut self, i: usize, j: usize, k: usize) { 
        self.insert_cob(i, j, k, |r0, r1| { 
            let c = CobComp::sdl(r0, r1);
            vec![c]
        })
    }

    fn insert_cob<F>(&mut self, i: usize, j: usize, k: usize, f: F) 
    where F: Fn((usize, usize), (usize, usize)) -> Vec<CobComp> { 
        let a = &self.mats[i][[k, j]]; // c[i][j] -> c[i+1][k]
        if a.is_zero() { 
            return
        }

        let u0 = &self.updates[i][j];
        let u1 = &self.updates[i+1][k];

        assert_eq!(u0.len(), 2);
        assert_eq!(u1.len(), 2);

        let i0 = u0[1].reindex(u0[0].index());
        let j0 = u0[1].index();
        let i1 = u1[1].reindex(u1[0].index());
        let j1 = u1[1].index();

        let r0 = (i0, j0);
        let r1 = (i1, j1);

        let a = std::mem::take(&mut self.mats[i][[k, j]]);
        let a = a.into_map_gens(|mut cob| { 
            cob.apply_update(&u0[0], End::Src);
            cob.apply_update(&u0[1], End::Src);
            cob.apply_update(&u1[0], End::Tgt);
            cob.apply_update(&u1[1], End::Tgt);

            for c in f(r0, r1).into_iter() { 
                cob.insert(c);
            }

            cob
        });

        self.mats[i][[k, j]] = a;
    }

    fn deloop(&mut self) { 
        for i in 0..self.objs.len() {
            let mut j = 0;
            while j < self.objs[i].len() { // number of objects changes by `deloop`.
                if let Some(r) = self.objs[i][j].find_loop() { 
                    self.deloop_at(i, j, r);
                } else { 
                    j += 1;
                }
            }
        }
    }

    fn deloop_at(&mut self, i: usize, j: usize, r: usize) { 
        let ci = &mut self.objs[i];

        let v0 = &mut ci[j];
        v0.deloop(r);

        let mut v1 = v0.clone();
        v0.append_label(KhAlgLabel::X);
        v1.append_label(KhAlgLabel::I);

        ci.insert(j+1, v1);

        if i > 0 { 
            self.cap_rows(i-1, j, r);
        }
        self.cup_cols(i, j, r);
    }

    fn cap_rows(&mut self, i: usize, j: usize, r: usize) { 
        let d = &mut self.mats[i];

        d.insert_zero_row(j+1);

        for k in 0..d.cols() { 
            let a0 = &d[[j, k]];
            if a0.is_zero() { 
                continue;
            }

            let a0 = std::mem::take(&mut d[[j, k]]);
            let a1 = a0.clone();

            d[[j,   k]] = a0.cap_off(r, Dot::None, End::Tgt);
            d[[j+1, k]] = a1.cap_off(r, Dot::Y,    End::Tgt);
        }
    }

    fn cup_cols(&mut self, i: usize, j: usize, r: usize) { 
        let d = &mut self.mats[i];

        d.insert_zero_col(j+1);
        
        for k in 0..d.rows() { 
            let a = &d[[k, j]];
            if a.is_zero() { 
                continue;
            }

            let a0 = std::mem::take(&mut d[[k, j]]);
            let a1 = a0.clone();

            d[[k, j]]   = a0.cap_off(r, Dot::X,    End::Src);
            d[[k, j+1]] = a1.cap_off(r, Dot::None, End::Src);
        }
    }
}

#[cfg(test)]
mod tests { 
    use yui_link::*;
    use yui_homology::*;
    use super::*;

    #[test]
    fn empty() { 
        let c = TngComplex::new();

        assert_eq!(c.len(), 1);
        assert_eq!(c.rank(0), 1);
    }

    #[test]
    fn single_x() { 
        let mut c = TngComplex::new();
        let x = Crossing::new(CrossingType::Xn, [0,1,2,3]);
        c.append(&x);

        assert_eq!(c.len(), 2);
        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 1);

        assert_eq!(c.tng(0, 0).ncomps(), 2);
        assert!(c.tng(0, 0).comp(0).is_arc());
        assert!(c.tng(0, 0).comp(1).is_arc());

        assert_eq!(c.tng(1, 0).ncomps(), 2);
        assert!(c.tng(1, 0).comp(0).is_arc());
        assert!(c.tng(1, 0).comp(1).is_arc());
    }

    #[test]
    fn two_x_disj() { 
        let mut c = TngComplex::new();
        let x0 = Crossing::new(CrossingType::Xn, [0,1,2,3]);
        let x1 = Crossing::new(CrossingType::Xn, [4,5,6,7]);

        c.append(&x0);
        c.append(&x1);

        assert_eq!(c.len(), 3);
        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 2);
        assert_eq!(c.rank(2), 1);

        assert_eq!(c.mat(0).shape(), (2, 1));
        assert_eq!(c.mat(1).shape(), (1, 2));
    }

    #[test]
    fn two_x() { 
        let mut c = TngComplex::new();
        let x0 = Crossing::new(CrossingType::Xn, [0,4,1,5]);
        let x1 = Crossing::new(CrossingType::Xn, [3,1,4,2]);

        c.append(&x0);
        c.append(&x1);

        assert_eq!(c.len(), 3);
        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 2);
        assert_eq!(c.rank(2), 1);

        assert_eq!(c.tng(0, 0).ncomps(), 2);
        assert_eq!(c.tng(1, 0).ncomps(), 2);
        assert_eq!(c.tng(1, 1).ncomps(), 2);
        assert_eq!(c.tng(2, 0).ncomps(), 3);
    }

    #[test]
    fn inv() { 
        let c = Cob::new(vec![
            CobComp::cyl(0, 1),
            CobComp::cyl(2, 3)
        ]);
        let f = Mor::from((c.clone(), -1));

        assert!(f.is_invertible());
        assert_eq!(f.inv(), Some(Mor::from(
            (Cob::new(vec![
                CobComp::cyl(1, 0),
                CobComp::cyl(3, 2)
            ]), -1)
        )));

        let f = Mor::from((c.clone(), 2));
        assert_eq!(f.is_invertible(), false);
        assert_eq!(f.inv(), None);
    }

    #[test]
    fn deloop_one() { 
        let mut c = TngComplex::new();
        let x0 = Crossing::new(CrossingType::Xp, [1,4,2,5]);
        let x1 = Crossing::new(CrossingType::Xn, [3,6,4,1]);

        c.append(&x0);
        c.append(&x1);

        assert_eq!(c.len(), 3);
        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 2);
        assert_eq!(c.rank(2), 1);

        assert!(c.tng(0, 0).find_loop().is_none());
        assert!(c.tng(1, 0).find_loop().is_some()); // loop here
        assert!(c.tng(1, 1).find_loop().is_none());
        assert!(c.tng(2, 0).find_loop().is_none());

        c.deloop();

        assert_eq!(c.len(), 3);
        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 3); // delooped here
        assert_eq!(c.rank(2), 1);

        assert!(c.tng(1, 0).find_loop().is_none()); // delooped
        assert!(c.tng(1, 1).find_loop().is_none()); // delooped

        assert_eq!(c.obj(1, 0).label, vec![KhAlgLabel::X]);
        assert_eq!(c.obj(1, 1).label, vec![KhAlgLabel::I]);
        
        c.describe()
    }

    #[test]
    fn trefoil_no_deloop() { 
        let mut c = TngComplex::new();
        let l = Link::trefoil();

        for x in l.data() {
            c.append(x);
        }

        assert_eq!(c.len(), 4);

        assert_eq!(c.rank(0), 1);
        assert_eq!(c.rank(1), 3);
        assert_eq!(c.rank(2), 3);
        assert_eq!(c.rank(3), 1);

        assert_eq!(c.mat(0).shape(), (3, 1));
        assert_eq!(c.mat(1).shape(), (3, 3));
        assert_eq!(c.mat(2).shape(), (1, 3));
        assert_eq!(c.mat(3).shape(), (0, 1));

        for i in 0..= 3 { 
            let ci = &c.objs[i];
            for j in 0..ci.len() { 
                let v = &ci[j];
                assert_eq!(v.state.weight(), i);
                assert!(v.tangle.is_closed());
            }
        }
    }

    #[test]
    fn trefoil_deloop() { 
        let mut c = TngComplex::new();
        let l = Link::trefoil();

        for x in l.data() {
            c.append(x);
        }

        c.deloop();

        assert_eq!(c.len(), 4);

        assert_eq!(c.rank(0), 8);
        assert_eq!(c.rank(1), 12);
        assert_eq!(c.rank(2), 6);
        assert_eq!(c.rank(3), 4);
    }

    #[test]
    fn trefoil_eval_homology() { 
        let mut c = TngComplex::new();
        let l = Link::trefoil();

        for x in l.data() {
            c.append(x);
        }

        c.deloop();

        let h = c.as_generic(&0, &0).homology(); // TODO should shift degree

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[0].is_free(), true);

        assert_eq!(h[1].rank(), 1);
        assert_eq!(h[1].tors(), &vec![2]);
        
        assert_eq!(h[2].is_zero(), true);
        assert_eq!(h[2].is_free(), true);
        
        assert_eq!(h[3].rank(), 2);
        assert_eq!(h[3].is_free(), true);
    }

    #[test]
    fn test() { 
        let mut c = TngComplex::new();
        let l = Link::trefoil();
        let mut step = 0;

        for x in l.data() {
            step += 1;
            c.append(x);

            println!("step: {step}");
            c.describe();
        }
    }
}
