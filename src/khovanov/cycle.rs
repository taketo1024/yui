use std::collections::HashSet;
use itertools::Itertools;
use crate::math::traits::{Ring, RingOps};
use crate::khovanov::algebra::KhAlgGen;
use crate::links::{links::{Component, State}, Link};
use super::complex::KhComplex;
use super::algebra::{KhAlgStr, KhEnhState, KhChain};

impl<R> KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn canon_cycle(&self) -> KhChain<R> {
        self._canon_cycle(true)
    }

    pub fn canon_cycles(&self) -> Vec<KhChain<R>> {
        vec![self._canon_cycle(true), self._canon_cycle(false)]
    }

    pub fn _canon_cycle(&self, positive: bool) -> KhChain<R> {
        assert!(self.structure().t().is_zero(), "Only supported for t = 0.");
        assert_eq!(self.link().components().len(), 1, "Only knots are supported.");

        let l = self.link();
        let s = l.ori_pres_state();
        let circles = l.resolved_by(&s).components();
        let colors = self.color_circles(l, &circles, positive);

        self.make_chain(&s, &colors)
    }

    fn color_circles(&self, l: &Link, circles: &Vec<Component>, positive: bool) -> Vec<Color> { 
        let n = circles.len();
        let mut colors = vec![Color::A; n];

        if n == 0 {
            return colors
        }

        let mut queue = vec![0];
        let mut remain: HashSet<_> = (1..n).collect();

        colors[0] = if positive { Color::A } else { Color::B };

        while !queue.is_empty() { 
            let i1 = queue.remove(0);
            let c1 = &circles[i1];

            let adjs = remain.iter().filter_map(|&i2| {
                let c2 = &circles[i2];
                if c1.is_adj(c2, l) { Some(i2) } else { None }
            }).collect_vec();
            
            for i2 in adjs {
                remain.remove(&i2);
                queue.push(i2);
                colors[i2] = colors[i1].other();
            };
        }

        colors
    }

    fn make_chain(&self, s: &State, colors: &Vec<Color>) -> KhChain<R> { 
        let str = self.structure();
        let mut z = KhChain::from((
            KhEnhState::new(s.clone(), vec![]),
            R::one()
        ));

        for c in colors {
            z *= c.as_tensor_factor(str);
        }

        z
    }
}

#[derive(PartialEq, Eq, Clone)]
enum Color { A, B }

impl Color { 
    fn other(&self) -> Self { 
        match self { 
            Color::A => Color::B,
            Color::B => Color::A
        }
    }

    fn as_tensor_factor<R>(&self, s: &KhAlgStr<R>) -> KhChain<R>
    where R: Ring, for<'x> &'x R: RingOps<R> {
        use KhAlgGen::{I, X};

        match self { 
            Color::A => KhChain::from((
                KhEnhState::new(State::empty(), vec![X]),
                R::one()
            )),
            Color::B => KhChain::from(vec![(
                KhEnhState::new(State::empty(), vec![X]),
                R::one()
            ), (
                KhEnhState::new(State::empty(), vec![I]),
                -s.h()
            )])
        }
    }
}

#[cfg(test)]
mod tests {
    use num_traits::Zero;

    use crate::{links::Link, khovanov::complex::KhComplex};
 
    #[test]
    fn trefoil() { 
        let l = Link::trefoil().mirror();
        let c = KhComplex::new_ht(l, 1, 0);
        let zs = c.canon_cycles();

        for z in zs { 
            let dz = c.differetiate(&z);
            assert_eq!(z.is_zero(), false);
            assert_eq!(dz.is_zero(), true);
        }
    }
 
    #[test]
    fn figure8() { 
        let l = Link::figure8();
        let c = KhComplex::new_ht(l, 1, 0);
        let zs = c.canon_cycles();

        for z in zs { 
            let dz = c.differetiate(&z);
            assert_eq!(z.is_zero(), false);
            assert_eq!(dz.is_zero(), true);
        }
    }
}