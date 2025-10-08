use yui::lc::{EitherGen, Lc};
use yui::poly::HPoly;
use yui::{Ring, RingOps, Sign};
use yui_homology::{ChainComplex, ChainMap};
use yui_link::Link;

use super::{KhAlgGen, KhAlgStr, KhChain, KhComplex, KhGen, KhLabel};

impl<R> KhAlgStr<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn sigma(&self, x: &KhAlgGen) -> Lc<KhAlgGen, R> {
        match x { 
            KhAlgGen::I => Lc::from((KhAlgGen::I, R::one())),
            KhAlgGen::X => Lc::from_iter([
                (KhAlgGen::X, -R::one()), 
                (KhAlgGen::I, self.h().clone())
            ])
        }
    }

    pub fn sigma_tensor(&self, x: &KhGen) -> KhChain<R> {
        let init = KhChain::from(
            KhGen::new(x.state, KhLabel::empty(), x.deg_shift)
        );

        x.label.iter().map(|x| self.sigma(&x)).fold(init, |res, next| { 
            res.combine(&next, |a, b| 
                KhGen::new(a.state, a.label + b, a.deg_shift)
            )
        })
    }

    pub fn sigma_chain(&self, z: &KhChain<R>) -> KhChain<R> {
        z.apply(|x| self.sigma_tensor(x))
    }
}

impl<R> KhAlgStr<HPoly<'H', R>>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn gr_sigma_chain(&self, z: &KhChain<HPoly<'H', R>>) -> KhChain<HPoly<'H', R>> {
        z.iter().flat_map(|(x, r)| {
            let deg = (r.deg() + x.label.iter().filter(|x| x.is_X()).count()) as i32;
            let e = HPoly::from_sign(Sign::from_parity(deg));
            let sx = self.sigma_tensor(x) * (e * r);
            sx.into_iter() 
        }).collect()
    }
}

impl<R> KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn sigma_complex(l: &Link, h: &R, t: &R) -> ChainComplex<EitherGen<KhGen, KhGen>, R> { 
        let c = KhComplex::new_no_simplify(&l, h, t, false);
        let (h0, h1) = c.h_range().into_inner();
        let str = c.str().clone();
        let inner = c.inner();

        let f = ChainMap::new(inner, inner, 0, move |_, z| z - str.sigma_chain(z));
        f.cone(inner, inner, h0 ..= h1 + 1, false)
    }
}

#[cfg(test)]
mod tests {
    use num_traits::Zero;
    use yui::poly::HPoly;
    use yui::{FF2};
    use yui_homology::{ChainComplexTrait, DisplaySeq};
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

        println!("{name}");
        c.homology().print_seq("i");
    }
}