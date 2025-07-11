use yui::lc::Lc;
use yui::poly::HPoly;
use yui::{Ring, RingOps, Sign};

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

    fn sigma_tensor(&self, x: &KhGen) -> KhChain<R> {
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
}

#[cfg(test)]
mod tests {
    use num_traits::Zero;
    use yui::poly::HPoly;
    use yui::{Ratio, FF2};
    use yui_homology::{ChainComplexTrait, ChainMap, DisplaySeq};
    use yui_link::Link;

    use super::*;

    #[test]
    fn test_sigma() {
        // type K = Ratio<i128>;
        type K = FF2;
        type R = HPoly<'H', K>;

        let h = R::variable();
        let t = R::zero();

        for i in 3..=7 { 
            for j in 1..=100 { 
                let name = format!("{}_{}", i, j);
                let Ok(l) = Link::load(&name) else { continue };
                let l = if l.writhe() >= 0 { l } else { l.mirror() };

                let c = KhComplex::new(&l, &h, &t, false);

                println!("{name}");
                c.homology().print_seq("i");
                
                let (h0, h1) = c.h_range().into_inner();
                let str = c.str().clone();
                let inner = c.inner();
        
                let f = ChainMap::new(inner, inner, 0, move |_, z| z - str.gr_sigma_chain(z));
                f.check_all(inner, inner);
        
                let cone = f.cone(inner, inner, h0 ..= h1 + 1, false);
                cone.check_d_all();
        
                let h = cone.homology();
                h.print_seq("i");
            }
        }
    }
}