use std::collections::HashSet;
use std::iter::zip;
use itertools::Itertools;
use yui_core::{Ring, RingOps};
use yui_link::{Link, LinkComp, State};

use crate::{KhAlgGen, KhComplex, KhLabel, KhEnhState, KhChain};

#[derive(PartialEq, Eq, Clone)]
pub enum Color { A, B }

impl Color { 
    pub fn is_a(&self) -> bool { 
        self == &Color::A
    }

    fn other(&self) -> Self { 
        match self { 
            Color::A => Color::B,
            Color::B => Color::A
        }
    }
}

pub trait ColoredSeifertCircles { 
    fn colored_seifert_circles(&self, positive: bool) -> Vec<(LinkComp, Color)>;
}

impl ColoredSeifertCircles for Link {
    fn colored_seifert_circles(&self, positive: bool) -> Vec<(LinkComp, Color)> {
        assert_eq!(self.components().len(), 1, "Only knots are supported.");

        let circles = self.seifert_circles();
        let n = circles.len();
    
        let mut colors = vec![Color::A; n];
        let mut queue = vec![];
        let mut remain: HashSet<_> = (0..n).collect();
    
        let e = self.first_edge().unwrap();
        let i = circles.iter().find_position(|c| c.edges().contains(e)).unwrap().0;
    
        queue.push(i);
        colors[i] = if positive { Color::A } else { Color::B };
    
        while !queue.is_empty() { 
            let i1 = queue.remove(0);
            let c1 = &circles[i1];
    
            let adjs = remain.iter().filter_map(|&i2| {
                let c2 = &circles[i2];
                if c1.is_adj(c2, self) { Some(i2) } else { None }
            }).collect_vec();
            
            for i2 in adjs {
                remain.remove(&i2);
                queue.push(i2);
                colors[i2] = colors[i1].other();
            };
        }
    
        assert!(queue.is_empty());
    
        zip(circles.into_iter(), colors.into_iter()).collect()
    }
}

impl<R> KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn canon_cycle(&self) -> KhChain<R> {
        self.canon_cycles().remove(0)
    }

    pub fn canon_cycles(&self) -> Vec<KhChain<R>> {
        let l = self.link();

        let str = self.structure();
        assert!(str.t().is_zero(), "Only supported for t = 0.");

        let a = R::zero(); // X_a = X
        let b = str.h();   // X_b = X - h
    
        let ori = if self.is_reduced() { 
            vec![true]
        } else { 
            vec![true, false]
        };

        ori.into_iter().map(|o| { 
            Self::make_canon_cycle(l, &a, b, o)
        }).collect()
    }

    fn make_canon_cycle(l: &Link, a: &R, b: &R, positive: bool) -> KhChain<R> {
        let s = l.ori_pres_state();
        let colors = l.colored_seifert_circles(positive);

        let mut z = KhChain::from(
            KhEnhState::new(s, KhLabel::empty())
        );

        let x_a = make_factor(a); // X - a
        let x_b = make_factor(b); // X - b
    
        for (_, c) in colors {
            z *= match c { 
                Color::A => &x_a,
                Color::B => &x_b
            }
        }
    
        z
    }
}

fn make_factor<R>(a: &R) -> KhChain<R> // a -> X - a
where R: Ring, for<'x> &'x R: RingOps<R> { 
    use KhAlgGen::{I, X};

    fn init(x: KhAlgGen) -> KhEnhState { 
        KhEnhState::new(
            State::empty(),
            KhLabel::from(x)
        )
    }

    KhChain::from_iter([
        (init(X), R::one()), 
        (init(I), -a)
    ])
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::Zero;
 
    #[test]
    fn trefoil() { 
        let l = Link::trefoil().mirror();
        let c = KhComplex::new(l, 1, 0, false);
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
        let c = KhComplex::new(l, 1, 0, false);
        let zs = c.canon_cycles();

        for z in zs { 
            let dz = c.differetiate(&z);
            assert_eq!(z.is_zero(), false);
            assert_eq!(dz.is_zero(), true);
        }
    }
}