use itertools::Either;
use yui::lc::{EitherGen, Lc};
use yui::poly::HPoly;
use yui::{EucRing, EucRingOps, Ring, RingOps, Sign};
use yui_homology::{ChainComplex, ChainMap, Grid2, Summand};
use yui_link::Link;

use crate::khi::KhIGen;
use crate::misc::make_gen_grid;

use super::{KhGen, KhAlg, KhChain, KhComplex, KhChainGen, KhTensor};

impl<R> KhAlg<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn sigma(&self, x: &KhGen) -> Lc<KhGen, R> {
        match x { 
            KhGen::I => Lc::from((KhGen::I, R::one())),
            KhGen::X => Lc::from_iter([
                (KhGen::X, -R::one()), 
                (KhGen::I, self.h().clone())
            ])
        }
    }

    pub fn sigma_tensor(&self, x: &KhChainGen) -> KhChain<R> {
        let init = KhChain::from(
            KhChainGen::new(x.state, KhTensor::empty(), x.deg_shift)
        );

        x.tensor.iter().map(|x| self.sigma(&x)).fold(init, |res, next| { 
            res.combine(&next, |a, b| 
                KhChainGen::new(a.state, a.tensor + b, a.deg_shift)
            )
        })
    }

    pub fn sigma_chain(&self, z: &KhChain<R>) -> KhChain<R> {
        z.apply(|x| self.sigma_tensor(x))
    }
}

impl<R> KhAlg<HPoly<'H', R>>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn gr_sigma_chain(&self, z: &KhChain<HPoly<'H', R>>) -> KhChain<HPoly<'H', R>> {
        z.iter().flat_map(|(x, r)| {
            let deg = (r.deg() + x.tensor.iter().filter(|x| x.is_X()).count()) as i32;
            let e = HPoly::from_sign(Sign::from_parity(deg));
            let sx = self.sigma_tensor(x) * (e * r);
            sx.into_iter() 
        }).collect()
    }
}

impl<R> KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn sigma_complex(l: &Link, h: &R, t: &R) -> ChainComplex<EitherGen<KhChainGen, KhChainGen>, R> { 
        let c = KhComplex::new_no_simplify(&l, h, t, false);
        let (h0, h1) = c.h_range().into_inner();
        let str = c.str().clone();
        let inner = c.inner();

        let f = ChainMap::new(inner, inner, 0, move |_, z| z - str.sigma_chain(z));
        f.cone(inner, inner, h0 ..= h1 + 1, false)
    }

    // TODO should unify KhIGen and EitherGen<KhGen, KhGen>.
    pub fn sigma_homology(l: &Link, h: &R, t: &R) -> Grid2<Summand<KhIGen, R>>
    where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
        let c = Self::sigma_complex(l, h, t);
        let h = c.homology();

        // convert EitherGen -> KhIGen
        let gens = h.map(|s| s.map_raw_gens(|x| { 
            match x.inner() {
                Either::Left(x)  => KhIGen::B(x.clone()),
                Either::Right(y) => KhIGen::Q(y.clone())
            }
        }));

        make_gen_grid(&gens)
    }
}

#[cfg(test)]
mod tests {
    use num_traits::Zero;
    use yui::poly::HPoly;
    use yui::{FF2};
    use yui_homology::{ChainComplexTrait, SummandTrait};
    use yui_link::Link;

    use super::*;

    #[test]
    fn test_sigma_cone() {
        type K = FF2;
        type R = HPoly<'H', K>;

        let h = R::variable();
        let t = R::zero();

        let name = "3_1";
        let l = Link::load(name).unwrap();
        let c = KhComplex::sigma_complex(&l, &h, &t);

        c.check_d_all();
    }


    #[test]
    fn test_sigma_homology() {
        type K = FF2;
        type R = HPoly<'H', K>;

        let h = R::variable();
        let t = R::zero();

        let name = "3_1";
        let l = Link::load(name).unwrap();
        let h = KhComplex::sigma_homology(&l, &h, &t);

        assert_eq!(h[(0, -1)].rank(), 1);
        assert_eq!(h[(1, -3)].rank(), 1);
        assert_eq!(h[(-2, -7)].tors().len(), 1);
        assert_eq!(h[(-2, -5)].tors().len(), 1);
        assert_eq!(h[(-1, -7)].tors().len(), 1);
        assert_eq!(h[(-1, -5)].tors().len(), 1);
        assert_eq!(h[( 1, -1)].tors().len(), 1);

        // use yui_homology::DisplayTable;
        // h.print_table("i", "j");
    }
}