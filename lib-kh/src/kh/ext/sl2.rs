
use std::collections::HashMap;
use std::fmt::Display;
use std::marker::PhantomData;
use std::ops::Index;

use itertools::Itertools;
use num_traits::Zero;
use yui_core::{AddMon, Field, FieldOps, RangeExt, Ring, RingOps, Sign};
use yui_homology::{GridTrait, SummandTrait, isize3};
use yui_link::Link;
use yui_matrix::MatTrait;
use yui_matrix::dense::snf::fnf;

use crate::ext::{Color, LinkExt};
use crate::kh::ext::cc::KhChainMap;
use crate::kh::internal::v1::cube::KhCube;
use crate::kh::{KhChain, KhChainExt, KhChainGen, KhComplex, KhHomology};

impl<R> KhComplex<R>
where
    R: Ring,
    for<'x> &'x R: RingOps<R>,
{
    pub fn sl2_map(&self, l: &Link) -> KhSl2Map<R> {
        assert!(l.is_knot());

        let cube = self.cube().clone();
        KhSl2Map::new(l, cube)
    }
}

pub struct KhSl2Map<R> where
    R: Ring,
    for<'x> &'x R: RingOps<R>
{ 
    path: Vec<(usize, Sign)>,
    cube: KhCube<R>
}

impl<R> KhSl2Map<R> where
    R: Ring,
    for<'x> &'x R: RingOps<R>
{ 
    pub fn new(l: &Link, cube: KhCube<R>) -> Self { 
        assert!(l.is_knot());

        let path = Self::make_path(l);
        Self { path, cube }
    }

    fn make_path(l: &Link) -> Vec<(usize, Sign)> {
        let mut res = vec![];
        let color = l.colored_seifert_circles(l.min_edge().unwrap());

        l.traverse_from((0, 0), |i, j| {
            if l.node(i).is_resolved() { return }

            let edge = l.node(i).edge(j);
            let c = color.iter().find(|(c, _)| c.contains(edge)).unwrap();
            let e = match c.1 {
                Color::A => Sign::Pos,
                Color::B => Sign::Neg
            };
            res.push((i, e));
        });

        res
    }

    pub fn h_deg(&self) -> isize{ -2 }
    pub fn q_deg(&self) -> isize{ -4 }

    fn apply_chi(&self, x: &KhChainGen, i: usize) -> KhChain<R> {
        if x.state[i].is_zero() {
            return KhChain::zero();
        }

        let t = x.state.edit(|s| s.set_0(i));
        self.cube.rev_d_to(x, &t, true)
    }

    fn apply_u(&self, x: &KhChainGen) -> KhChain<R> { 
        let n = self.path.len();
        if n < 2 { 
            return KhChain::zero();
        }

        let indices = (0 .. n - 1).flat_map(|j1| {
            (j1 + 1 .. n).map(move |j2| (j1, j2) )
        });

        KhChain::sum(indices.map(|(j1, j2)| { 
            let (i1, e1) = self.path[j1];
            let (i2, e2) = self.path[j2];
            KhChain::from(*x).apply(|x|
                self.apply_chi(x, i1) * R::from_sign(e1)
            ).apply(|x|
                self.apply_chi(x, i2) * R::from_sign(e2)
            )
        }))
    }

    pub fn apply(&self, z: &KhChain<R>) -> KhChain<R> { 
        z.apply(|x| self.apply_u(x))
    }

    pub fn into_chain_map(self) -> KhChainMap<R> { 
        KhChainMap::new(self.h_deg(), move |_, z| {
            self.apply(z)
        })
    }

    pub fn string_decomp(&self, kh: &KhHomology<R>) -> StringDecomp<R>
    where R: Field, for<'x> &'x R: FieldOps<R> {
        use yui_matrix::sparse::SpMat;
        
        assert!(self.cube.str().h().is_zero());
        assert!(self.cube.str().t().is_zero());

        let h_range = kh.h_range().mv(0, -self.h_deg());

        let n = kh.support().map(|&i| kh[i].rank()).sum();
        let gens = h_range.clone().flat_map(|i| kh[i].gens()).collect_vec();

        let blocks = h_range.map(|i|
            kh[i].make_matrix(&kh[i + self.h_deg()], |z| self.apply(z))
        );
        let e_mat = SpMat::block_diag(blocks);

        assert_eq!(gens.len(), n);
        assert_eq!(e_mat.shape(), (n, n));

        let flags = [false, true, false, false]; // only need pinv
        let fnf = fnf(&e_mat.clone().into_dense(), flags);
        let (res, [_, pinv, ..]) = fnf.destruct();
        let pinv = pinv.unwrap().into_sparse();

        let data = (0..n).flat_map(|i| { 
            let ord = res[(i, i)].lead_deg(); // extract torsion order l from x^l.
            if ord == 0 { return None; }

            let v = pinv.col_vec(i);

            // println!("{i}) order: {ord}\n{:?}", v.to_dense());

            let indices = v.iter_nz().filter_map(|(j, r)| 
                if r.is_const() { 
                    Some(j) 
                } else { 
                    None 
                }
            ).collect_vec();

            assert!(!indices.is_empty());
            assert!(indices.iter().map(|&j| gens[j].q_deg()).all_equal());

            let j = *indices.first().unwrap();
            let z = &gens[j];
            let t = z.h_deg();
            let q = -z.q_deg(); // TODO: must modify q_deg later.
            let d = 2 * t + q;

            Some(isize3(d, q, ord as isize))
        }).counts();

        StringDecomp::new(data)
    }
}

#[derive(Debug, Clone)]
pub struct StringDecomp<R> where
    R: Ring,
    for<'x> &'x R: RingOps<R>
{ 
    data: HashMap<isize3, usize>, // {(i, j, l) : n} => t^i q^j (R[x]/(x^l))^n
    _phantom: PhantomData<R>
}

impl<R> StringDecomp<R> where
    R: Ring,
    for<'x> &'x R: RingOps<R>
{
    fn new(data: HashMap<isize3, usize>) -> Self { 
        Self { data, _phantom: PhantomData }
    } 

    pub fn len(&self) -> usize {
        self.data.len()
    }
}

impl<R> Index<(isize, isize, isize)> for StringDecomp<R>
where
    R: Ring,
    for<'x> &'x R: RingOps<R>,
{
    type Output = usize;

    fn index(&self, index: (isize, isize, isize)) -> &Self::Output {
        self.data.get(&index.into()).unwrap_or(&0)
    }
}

impl<R> From<HashMap<isize3, usize>> for StringDecomp<R>
where
    R: Ring,
    for<'x> &'x R: RingOps<R>
{
    fn from(data: HashMap<isize3, usize>) -> Self {
        Self::new(data)
    }
}

impl<R> FromIterator<(isize3, usize)> for StringDecomp<R>
where
    R: Ring,
    for<'x> &'x R: RingOps<R>
{
    fn from_iter<T>(iter: T) -> Self
    where T: IntoIterator<Item = (isize3, usize)> {
        let data = HashMap::from_iter(iter);
        Self::new(data)
    }
}

impl<R> Display for StringDecomp<R>
where
    R: Ring,
    for<'x> &'x R: RingOps<R>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use yui_core::util::format;
        use yui_core::poly::Var2;
        type V = Var2<'δ', 'q', isize>;

        let list = self.data.iter().sorted_by_key(|(i, _)| *i);
        let str = format::lc(list.map(|(&i, n)| {
            let (d, q, l) = i.into();
            let v = V::from((d, q));
            let t = if v.multi_deg() == (0, 0) { 
                format!("e({l})")
            } else { 
                format!("{}e({l})", v)
            };
            (t, n)
        }));
        write!(f, "{str}")
    }
}

#[cfg(test)]
mod tests {
    use yui_core::num::Ratio;
    use yui_link::State;
    
    #[allow(unused)]
    use yui_homology::DisplayTable;

    use crate::kh::KhChainExt;

    use super::*;

    #[test]
    fn test_sl2map_unknot() {
        let l = Link::unknot();
        let cube = KhCube::new(&l, &0, &0, None, (0, 0));
        let map = KhSl2Map::new(&l, cube);

        assert_eq!(map.path.len(), 0);
    } 

    #[test]
    fn test_sl2map_2twist_unknot() {
        let l = Link::from_pd_code([[1,1,2,4],[3,3,4,2]]);
        let cube = KhCube::new(&l, &0, &0, None, (0, 0));
        let map = KhSl2Map::new(&l, cube);

        assert_eq!(map.path.len(), 4);
        assert_eq!(map.path, vec![(0, Sign::Pos), (1, Sign::Neg), (1, Sign::Pos), (0, Sign::Neg)])
    } 

    #[test]
    fn test_sl2map_trefoil() {
        let l = Link::trefoil();
        let cube = KhCube::new(&l, &0, &0, None, (0, 0));
        let map = KhSl2Map::new(&l, cube);

        assert_eq!(map.path.len(), 6);
        assert_eq!(map.path, vec![(0, Sign::Pos), (2, Sign::Neg), (1, Sign::Pos), (0, Sign::Neg), (2, Sign::Pos), (1, Sign::Neg)])
    } 

    #[test]
    fn test_u_unknot() {
        let l = Link::unknot();
        let cube = KhCube::new(&l, &0, &0, None, (0, 0));
        let map = KhSl2Map::new(&l, cube.clone());

        let v = State::empty();
        let x = cube.vertex(&v).generators()[0];
        assert_eq!(map.apply_u(x), KhChain::zero());
    } 

    #[test]
    fn test_u_2twist_unknot() {
        let l = Link::from_pd_code([[1,1,2,4],[3,3,4,2]]);
        let cube = KhCube::new(&l, &0, &0, None, (0, 0));
        let map = KhSl2Map::new(&l, cube.clone());

        let v = State::from([0,0]);
        let x = cube.vertex(&v).generators()[0]; // 111
        assert_eq!(map.apply_u(x), KhChain::zero());

        let v = State::from([1,1]);
        let x = cube.vertex(&v).generators()[0]; // 1
        assert_eq!(map.apply_u(x), KhChain::zero());
    } 

    #[test]
    fn test_u_trefoil() {
        let l = Link::trefoil();
        let cube = KhCube::new(&l, &0, &0, None, (0, 0));
        let map = KhSl2Map::new(&l, cube.clone());

        let v = State::from([1,1,1]);
        let z = cube.vertex(&v).generators()[0]; // (11)₁₁₁
        let w = map.apply_u(z);

        assert_ne!(w, KhChain::zero());

        assert_eq!(z.h_deg(), 3);
        assert_eq!(z.q_deg(), 5);
        assert_eq!(w.h_deg(), 1);
        assert_eq!(w.q_deg(), 1);
    } 

    #[test]
    fn test_ch_map_unknot() {
        let l = Link::unknot();
        let c = KhComplex::new_no_simplify(&l, &0, &0, false);
        let e = c.sl2_map(&l).into_chain_map();

        assert_eq!(e.deg(), -2);
        e.check_all(c.inner(), c.inner());
    }

    #[test]
    fn test_ch_map_2twist_unknot() {
        let l = Link::from_pd_code([[1,1,2,4],[3,3,4,2]]);
        let c = KhComplex::new_no_simplify(&l, &0, &0, false);
        let e = c.sl2_map(&l).into_chain_map();

        assert_eq!(e.deg(), -2);
        e.check_all(c.inner(), c.inner());
    } 

    #[test]
    fn test_ch_map_trefoil() {
        let l = Link::trefoil();
        let c = KhComplex::new_no_simplify(&l, &0, &0, false);
        let e = c.sl2_map(&l).into_chain_map();

        assert_eq!(e.deg(), -2);
        e.check_all(c.inner(), c.inner());

        let z = c[0].gen(0); // (11)₁₁₁
        let w = e.apply(0, &z);

        assert_ne!(w, KhChain::zero());
        
        assert_eq!(z.h_deg(), 0);
        assert_eq!(z.q_deg(), -1);
        assert_eq!(w.h_deg(), -2);
        assert_eq!(w.q_deg(), -5);
    }

    #[test]
    fn test_ch_map_6_1() {
        let l = Link::load("6_1").unwrap();
        let c = KhComplex::new_no_simplify(&l, &0, &0, false);
        let e = c.sl2_map(&l).into_chain_map();

        e.check_all(c.inner(), c.inner());
    }

    #[test]
    fn test_ch_map_6_2() {
        let l = Link::load("6_2").unwrap();
        let c = KhComplex::new_no_simplify(&l, &0, &0, false);
        let e = c.sl2_map(&l).into_chain_map();

        e.check_all(c.inner(), c.inner());
    }

    #[test]
    fn test_ch_map_6_3() {
        let l = Link::load("6_3").unwrap();
        let c = KhComplex::new_no_simplify(&l, &0, &0, false);
        let e = c.sl2_map(&l).into_chain_map();

        e.check_all(c.inner(), c.inner());
    }

    type QQ = Ratio<i64>;

    #[test]
    fn test_string_decomp_3_1() { 
        let l = Link::load("3_1").unwrap();
        let c = KhComplex::new_no_simplify(&l, &QQ::zero(), &QQ::zero(), true);
        let e = c.sl2_map(&l);
        let h = c.homology();

        let e_str = e.string_decomp(&h);

        assert_eq!(e_str.len(), 2);
        assert_eq!(e_str[(2, 8, 1)], 1);
        assert_eq!(e_str[(2, 2, 2)], 1);

        assert_eq!(e_str.to_string(), String::from("δ²q²e(2) + δ²q⁸e(1)"))
    }

    #[test]
    fn test_string_decomp_unred_3_1() { 
        let l = Link::load("3_1").unwrap();
        let c = KhComplex::new_no_simplify(&l, &QQ::zero(), &QQ::zero(), false);
        let e = c.sl2_map(&l);
        let h = c.homology();

        let e_str = e.string_decomp(&h);
        assert_eq!(e_str.len(), 3);
        assert_eq!(e_str[(1, 1, 2)], 1);
        assert_eq!(e_str[(3, 3, 1)], 1);
        assert_eq!(e_str[(3, 9, 1)], 1);

        assert_eq!(e_str.to_string(), String::from("δqe(2) + δ³q³e(1) + δ³q⁹e(1)"))
    }
}