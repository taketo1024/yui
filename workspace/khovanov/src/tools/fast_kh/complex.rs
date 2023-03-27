#![allow(dead_code)] // TODO remove 

use std::fmt::Display;

use itertools::Itertools;
use yui_lin_comb::LinComb;
use yui_matrix::sparse::MatType;
use yui_polynomial::Poly2;
use yui_matrix::dense::Mat;
use yui_link::{Crossing, Resolution, State, Component};

use crate::KhAlgLabel;

use super::cob::{Cob, Dot};
use super::tng::Tng;

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

    fn append_arc(&mut self, a: Component) {
        self.tangle.append_arc(a);
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
        self.objs.len() - 1
    }

    pub fn rank(&self, i: usize) -> usize { 
        if i < self.objs.len() as usize { 
            self.objs[i as usize].len()
        } else { 
            0
        }
    }

    pub fn obj(&self, i: usize, j: usize) -> &Obj { 
        &self.objs[i][j]
    }

    pub fn mat(&self, i: usize) -> &Mat<Mor> { 
        &self.mats[i]
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

    pub fn append_x(&mut self, x: &Crossing) {
        let n = self.objs.len();
        let ranks = (0 .. n+1).map(|i| self.rank(i)).collect_vec();

        self.extend_objs();
        self.extend_mats();

        // update tangles
        for i in 0..n { 
            let r = ranks[i];
            let s = ranks[i+1];
            for j in 0..r { 
                self.append_x_res(i,   j,   x, Resolution::Res0);
                self.append_x_res(i+1, s+j, x, Resolution::Res1);
            }
        }

        // update cobs
        for i in 0..n { 
            let r = ranks[i];
            let s = ranks[i+1];
            for j in 0..r { 
                // TODO: update blocks
            }
        }
    }

    fn extend_objs(&mut self) {
        let n = self.objs.len();
        let mut clone = self.objs.clone();

        self.objs.push(vec![]); // (n+1)

        for i in 0..n { 
            self.objs[i+1].append(&mut clone[i]);
        }
    }

    fn extend_mats(&mut self) { 
        let n = self.mats.len();
        let mut mats = vec![];

        for i in 0..n { 
            let d0 = &self.mats[i];
            let (l, k) = d0.shape();
            
            let d1 = if i > 0 {  
                -&self.mats[i-1]
            } else { 
                Mat::zero((k, 0))
            };

            let o = Mat::zero((l, d1.cols()));
            let f = Mat::diag((k, k), (0..k).map(|j| {
                let v = &self.objs[i][j];
                let f = Cob::id_for(&v.tangle);
                Mor::from(f)
            }));
            let d = Mat::combine_blocks([
                &d0, &o, 
                &f,  &d1
            ]);

            mats.push(d)
        }

        let d_n = Mat::zero((0, self.mats[n-1].cols()));
        mats.push(d_n);

        self.mats = mats
    }

    fn append_x_res(&mut self, i: usize, j: usize, x: &Crossing, r: Resolution) {
        let v = &mut self.objs[i][j];
        v.append_bit(r);

        let (r0, r1) = x.res_arcs(Resolution::Res0);
        self.append_arc(i, j, r0);
        self.append_arc(i, j, r1);
    }

    fn append_arc(&mut self, i: usize, j: usize, arc: Component) {
        let v = &mut self.objs[i][j];
        v.append_arc(arc);

        // TODO update cob
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

#[cfg(test)]
mod tests { 
    use yui_link::Link;
    use super::*;

    #[test]
    fn append_link() { 
        let mut c = TngComplex::new();
        let l = Link::trefoil();

        for x in l.data() {
            c.append_x(x);
        }

        assert_eq!(c.len(), 3);

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
}
