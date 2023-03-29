#![allow(unused)] // TODO remove 

use std::fmt::Display;

use itertools::Itertools;
use num_traits::Zero;
use cartesian::cartesian;
use yui_lin_comb::LinComb;
use yui_matrix::sparse::MatType;
use yui_polynomial::Poly2;
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

type R = Poly2<'H', 'T', i64>; // R = Z[H, T]
type Mor = LinComb<Cob, R>;    // Mor = Linear combinations of dotted-cobordisms over Z[H].

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
        let objs = std::mem::take(&mut self.objs);
        let mats = std::mem::take(&mut self.mats);

        let mut builder = TngComplexBuilder::new(objs, mats);
        builder.append(x);

        (self.objs, self.mats) = builder.result();
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

        let i0 = u0[1].apply(u0[0].index());
        let j0 = u0[1].index();
        let i1 = u1[1].apply(u1[0].index());
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
}

#[cfg(test)]
mod tests { 
    use yui_link::*;
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

        c.describe();

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
