
use num_traits::Zero;
use yui_core::{AddMon, Ring, RingOps, Sign};
use yui_link::Link;

use crate::kh::ext::cc::KhChainMap;
use crate::kh::internal::v1::cube::KhCube;
use crate::kh::{KhChain, KhChainGen, KhComplex};

impl<R> KhComplex<R>
where
    R: Ring,
    for<'x> &'x R: RingOps<R>,
{
    pub fn e_map(&self, l: &Link) -> KhChainMap<R> {
        assert!(l.is_knot());

        let cube = self.cube().clone();
        let map = KhSl2Map::new(l, cube);

        KhChainMap::new(self.inner(), self.inner(), -2, move |_, z| {
            z.apply(|x| map.apply_u(x))
        })
    }
}

struct KhSl2Map<R> where
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
    fn new(l: &Link, cube: KhCube<R>) -> Self { 
        let path = Self::make_path(l);
        Self { path, cube }
    }

    fn make_path(l: &Link) -> Vec<(usize, Sign)> {
        let mut res = vec![];
        let mut data = vec![0; l.n_nodes()];

        l.traverse_from((0, 0), |i, j| {
            if l.node(i).is_resolved() { return }
            data[i] += j;
        });

        l.traverse_from((0, 0), |i, j| {
            if l.node(i).is_resolved() { return }
            let e = match (data[i], j) { 
                (1, 0) | (3, 1) | (5, 2) | (3, 3) => Sign::Pos, // Left
                (1, 1) | (3, 2) | (5, 3) | (3, 0) => Sign::Neg, // Right
                _ => panic!("Invalid knot data.")
            };
            res.push((i, e));
        });

        res
    }

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
}

#[cfg(test)]
mod tests {
    use yui_link::State;

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
        assert_eq!(map.path, vec![(0, Sign::Neg), (1, Sign::Pos), (1, Sign::Neg), (0, Sign::Pos)])
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
    fn test_unknot() {
        let l = Link::unknot();
        let c = KhComplex::new_no_simplify(&l, &0, &0, false);
        let e = c.e_map(&l);

        e.check_all(c.inner(), c.inner());
    }

    #[test]
    fn test_trefoil() {
        let l = Link::trefoil();
        let c = KhComplex::new_no_simplify(&l, &0, &0, false);
        let e = c.e_map(&l);

        e.check_all(c.inner(), c.inner());
    }
}
