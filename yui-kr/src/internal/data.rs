use std::collections::{HashSet, HashMap};
use std::ops::RangeInclusive;
use std::sync::Arc;
use cartesian::{cartesian, TuplePrepend};
use itertools::{Itertools, izip};
use num_integer::Integer;
use num_traits::One;

use yui::{Ring, RingOps, Sign, AddMon};
use yui_homology::{isize3, GridIter};
use yui_link::util::seifert_graph;
use yui_link::{Link, CrossingType, Crossing, Edge};
use yui::bitseq::{BitSeq, Bit};

use super::base::KRPoly;
use super::hor_excl::KRHorExcl;

/*
 *    a   b
 *     \ /
 *      X
 *     / \
 *    c   d
 * 
 * x_ac = x_a - x_c = -(x_b - x_d),
 * x_bc = x_b - x_c = -(x_a - x_d)
 */

#[derive(Debug, Clone)]
pub struct KRCrossData<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub x_ac: KRPoly<R>,
    pub x_bc: KRPoly<R>
}

#[derive(Clone)]
pub struct KRCubeData<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    n_cross: usize,
    writhe: isize, 
    n_seif: isize,
    root_grad: isize3,
    x_signs: Vec<Sign>,
    x_polys: Vec<KRCrossData<R>>,
    verts: HashMap<usize, Vec<BitSeq>>,
    excl: HashMap<BitSeq, Arc<KRHorExcl<R>>>
}

impl<R> KRCubeData<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new(link: &Link) -> Self {
        Self::new_excl(link, 2)
    }

    pub(crate) fn new_excl(link: &Link, excl_level: usize) -> Self {
        assert!(excl_level <= 2);
        let mut data = Self::new_no_excl(link);
        data.perform_excl(excl_level);
        data
    }

    pub(crate) fn new_no_excl(link: &Link) -> Self { 
        assert!(link.is_knot(), "only knots are supported.");

        let n = link.crossing_num();
        let w = link.writhe() as isize;
        let s = link.seifert_circles().len() as isize;
        let x_signs = link.crossing_signs();
        let x_polys = Self::collect_x_polys(link);
        let root_grad = Self::grad(n, w, s, &x_signs, BitSeq::zeros(n), BitSeq::zeros(n));
        let verts = BitSeq::generate(n).into_group_map_by(|b| b.weight());
        let excl = HashMap::new();

        Self {
            n_cross: n,
            writhe: w,
            n_seif: s,
            root_grad,
            x_signs,
            x_polys,
            verts,
            excl // empty
        }
    }

    fn perform_excl(&mut self, excl_level: usize) {
        let n = self.n_cross;
        self.excl = BitSeq::generate(n).map(|v| {
            let excl = Arc::new( KRHorExcl::from(&self, v, excl_level) );
            (v, excl)
        }).collect();
    }

    pub fn dim(&self) -> usize { 
        self.n_cross
    }

    pub fn root_grad(&self) -> isize3 { 
        self.root_grad
    }

    pub fn x_sign(&self, i: usize) -> Sign { 
        self.x_signs[i]
    }

    pub fn x_poly(&self, i: usize) -> &KRCrossData<R> {
        &self.x_polys[i]
    }

    pub fn excl(&self, v_coods: BitSeq) -> Arc<KRHorExcl<R>> { 
        Arc::clone( &self.excl[&v_coods] )
    }

    pub fn verts_of_weight(&self, k: usize) -> &[BitSeq] {
        if let Some(vs) = self.verts.get(&k) { 
            &vs
        } else {
            &[]
        }
    }

    fn grad(n: usize, w: isize, s: isize, x_signs: &Vec<Sign>, h_coords: BitSeq, v_coords: BitSeq) -> isize3 { 
        assert_eq!(h_coords.len(), n);
        assert_eq!(v_coords.len(), n);

        let init = Self::l_grad_shift(w, s);
        let trip = izip!(
            x_signs.iter(), 
            h_coords.iter(), 
            v_coords.iter()
        );

        trip.fold(init, |grad, (&e, h, v)| { 
            grad + Self::x_grad_shift(e, h, v)
        })
    }

    fn l_grad_shift(w: isize, s: isize) -> isize3 { 
        isize3(-w + s - 1, w + s - 1, w - s + 1)
    }

    fn x_grad_shift(x_sign: Sign, h: Bit, v: Bit) -> isize3 { 
        use Bit::{Bit0, Bit1};
        if x_sign.is_positive() {
            match (h, v) {
                (Bit0, Bit0) => isize3(2, -2, -2),
                (Bit1, Bit0) => isize3(0,  0, -2),
                (Bit0, Bit1) => isize3(0, -2,  0),
                (Bit1, Bit1) => isize3(0,  0,  0)
            }
        } else { 
            match (h, v) {
                (Bit0, Bit0) => isize3( 0, -2, 0),
                (Bit1, Bit0) => isize3( 0,  0, 0),
                (Bit0, Bit1) => isize3( 0, -2, 2),
                (Bit1, Bit1) => isize3(-2,  0, 2)
            }
        }
    }

    pub fn grad_at(&self, h_coords: BitSeq, v_coords: BitSeq) -> isize3 { 
        Self::grad(self.n_cross, self.writhe, self.n_seif, &self.x_signs, h_coords, v_coords)
    }

    // Morton bound (only for knots)
    pub fn i_range(&self) -> RangeInclusive<isize> { 
        let n = self.n_cross as isize;
        let s = self.n_seif;
        -n + s - 1 ..= n - s + 1
    }

    // MFW inequality
    pub fn j_range(&self) -> RangeInclusive<isize> { 
        let w = self.writhe;
        let s = self.n_seif;
        w - s + 1 ..= w + s - 1
    }

    pub fn k_range(&self) -> RangeInclusive<isize> { 
        self.i_range()
    }

    pub fn q_range(&self) -> RangeInclusive<isize> { 
        range(self.support().filter_map(|idx| {
            self.to_inner_grad(idx).map(|g| g.0)
        }))
    }

    pub fn hv_range(&self, q: isize) -> (RangeInclusive<isize>, RangeInclusive<isize>) { 
        let hv = self.support().filter_map(|idx| {
            let inner = self.to_inner_grad(idx).unwrap();
            (inner.0 == q).then(|| (inner.1, inner.2))
        });
        range2(hv)
    }

    pub fn support(&self) -> GridIter<isize3> {
        let i_range = self.i_range().step_by(2);
        let j_range = self.j_range().step_by(2);
        let k_range = self.k_range().step_by(2);
        let range = cartesian!(i_range, j_range.clone(), k_range.clone());
        
        range.map(isize3::from).filter(|&idx| 
            self.is_supported(idx)
        ).collect_vec().into_iter()
    }

    pub fn is_supported(&self, grad: isize3) -> bool { 
        let isize3(i, j, k) = grad;
        let isize3(i0, j0, k0) = self.root_grad;

        let i_range = self.i_range();
        let j_range = self.j_range();
        let k_range = self.k_range();

        (i - i0).is_even() 
            && (j - j0).is_even() 
            && (k - k0).is_even()
            && i_range.contains(&i)
            && j_range.contains(&j)
            && k_range.contains(&k)
            && k_range.contains(&(k + 2 * i))
    }

    pub fn to_inner_grad(&self, grad: isize3) -> Option<isize3> { 
        let isize3(i, j, k) = grad;
        if i > 0 { 
            return self.to_inner_grad(isize3(-i, j, k + 2 * i))
        }

        let isize3(i0, j0, k0) = self.root_grad;
        
        if (i - i0).is_even() && (j - j0).is_even() && (k - k0).is_even() {
            let q = ((i - i0) - (j - j0)) / 2;
            let h = (j - j0) / 2;
            let v = (k - k0) / 2;
            Some(isize3(q, h, v))
        } else { 
            None
        }
    }

    pub fn to_outer_grad(&self, grad: isize3) -> isize3 { 
        let isize3(q, h, v) = grad;
        let isize3(i0, j0, k0) = self.root_grad;
        let i = 2 * q + 2 * h + i0;
        let j = 2 * h + j0;
        let k = 2 * v + k0;
        isize3(i, j, k)
    }

    pub fn hor_edge_poly(&self, v_coords: BitSeq, i: usize) -> KRPoly<R> {
        use Bit::{Bit0, Bit1};

        let sign = self.x_sign(i);
        let p = self.x_poly(i);
        let v = v_coords[i];

        match (sign.is_positive(), v) {
            (true, Bit0) | (false, Bit1) => &p.x_ac * &p.x_bc,
            (true, Bit1) | (false, Bit0) => p.x_ac.clone()
        }
    }

    pub fn ver_edge_poly(&self, h_coords: BitSeq, i: usize) -> KRPoly<R> {
        use Bit::{Bit0, Bit1};

        let sign = self.x_sign(i);
        let p = self.x_poly(i);
        let h = h_coords[i];

        match (sign.is_positive(), h) {
            (true, Bit0) | (false, Bit1) => p.x_bc.clone(),
            (true, Bit1) | (false, Bit0) => KRPoly::one()
        }
    }

    fn collect_x_polys(link: &Link) -> Vec<KRCrossData<R>> {
        let n = link.crossing_num();
        let signs = link.crossing_signs();
        let resolved = Self::resolve(link);

        assert_eq!(resolved.components().len(), 1); 

        let (path, mons) = { 
            let mut path = vec![];
            let mut mons = vec![];

            resolved.traverse_edges((0, 0), |i, j| { 
                let x = &resolved.data()[i];
                let e = x.edge(j);
                let x_i = KRPoly::variable(i);
                let p_i = match (Self::is_vert(x, signs[i]), j) {
                    (true, 0) | (false, 3) =>  x_i,
                    (true, 1) | (false, 0) => -x_i,
                    _ => panic!()
                };

                path.push(e);
                mons.push(p_i);
            });

            assert_eq!(path.len(), 2 * n + 1);
            assert_eq!(path.first(), path.last());

            path.pop();
            mons.pop();

            (path, mons)
        };

        let traverse = |from: Edge, to: Edge| -> KRPoly<R> { 
            let m = path.len(); // == 2 * n
            let i = path.iter().find_position(|&&e| e == from).unwrap().0;
            let j = path.iter().find_position(|&&e| e == to  ).unwrap().0;
            let l = if i < j { 
                j - i
            } else { 
                m + j - i
            };
            KRPoly::sum(
                (i .. i + l).map(|k| &mons[k % m])
            )
        };

        (0 .. n).map(|i| {
            let x = &resolved.data()[i];
            let (a,b,c) = if Self::is_vert(x, signs[i]) { 
                (3,2,0)
            } else { 
                (2,1,3)
            };

            let x_ac = traverse(x.edge(c), x.edge(a));
            let x_bc = traverse(x.edge(c), x.edge(b));
    
            KRCrossData { x_ac, x_bc }
        }).collect()
    }

    fn resolve(link: &Link) -> Link {
        use petgraph::algo::min_spanning_tree;
        use petgraph::data::Element::Edge as GraphEdge;

        let g = seifert_graph(link);
        let t = min_spanning_tree(&g).filter_map(|e| 
            match e { 
                GraphEdge { weight: x, ..} => Some(x),
                _ => None
            }
        ).collect::<HashSet<_>>();

        let n = link.crossing_num();
        let s0 = link.ori_pres_state();
        let mut data = link.data().clone();

        for i in 0..n { 
            if t.contains(&i) { 
                continue
            }
            data[i].resolve(s0[i])
        }

        Link::new(data)
    }

    // is_vert: 
    //          
    //  a(3) b(2)   c(3) a(2)
    //    ↑\ /↑        \ /→  
    //      X           X    
    //     / \         / \→   
    //  c(0) d(1)   d(0) b(1)
    //
    //    true        false

    fn is_vert(x: &Crossing, sign: Sign) -> bool { 
        use CrossingType::*;
        match (x.ctype(), sign.is_positive()) {
            (X, false) | (Xm, true)  | (V, _) => true,
            (X, true)  | (Xm, false) | (H, _) => false
        }
    }
}

pub fn range<I>(vals: I) -> RangeInclusive<isize>
where I: IntoIterator<Item = isize> { 
    let set = vals.into_iter().collect::<HashSet<_>>();
    let i0 = set.iter().min().cloned().unwrap_or(0);
    let i1 = set.iter().max().cloned().unwrap_or(0);
    i0..=i1
}

pub fn range2<I>(vals: I) -> (RangeInclusive<isize>, RangeInclusive<isize>)
where I: IntoIterator<Item = (isize, isize)> { 
    let set = vals.into_iter().collect::<HashSet<_>>();
    let r0 = range(set.iter().map(|i| i.0));
    let r1 = range(set.iter().map(|i| i.1));
    (r0, r1)
}

#[cfg(test)]
mod tests {
    use yui::Ratio;
    use yui_link::Link;

    use super::*;

    type R = Ratio<i64>;
    type P = KRPoly<R>;

    #[test]
    fn grad_shift() { 
        let l = Link::trefoil(); // (w, s) = (-3, 2)
        let w = l.writhe() as isize;
        let s = l.seifert_circles().len() as isize;
        assert_eq!(KRCubeData::<i32>::l_grad_shift(w, s), isize3(4,-2,-4));
    }

    #[test]
    fn grad_at() { 
        let l = Link::trefoil();
        let data = KRCubeData::<R>::new_excl(&l, 0);
        let grad0 = data.root_grad();
        assert_eq!(grad0, isize3(4,-8,-4));

        let h = BitSeq::from([1,0,0]);
        let v = BitSeq::from([0,1,0]);
        let grad1 = data.grad_at(h, v);
        assert_eq!(grad1, isize3(4,-6,-2));
    }

    #[test]
    fn resolve() { 
        let l = Link::hopf_link();
        let u = KRCubeData::<R>::resolve(&l);
        assert_eq!(u.components().len(), 1);
    }

    #[test]
    fn collect_x_polys() { 
        let l = Link::from_pd_code([[1,4,2,5],[5,2,6,3],[3,6,4,1]]); // trefoil
        let g = KRCubeData::<R>::collect_x_polys(&l);
        let x = (0..3).map(P::variable).collect_vec();

        assert_eq!(g.len(), 3);
        assert_eq!(g[0].x_ac, -&x[1] - &x[2]);
        assert_eq!(g[0].x_bc, x[0]);
        assert_eq!(g[1].x_ac, x[1]);
        assert_eq!(g[1].x_bc, &x[0] + &x[2]);
        assert_eq!(g[2].x_ac, x[2]);
        assert_eq!(g[2].x_bc, &x[0] - &x[1]);
    }
}