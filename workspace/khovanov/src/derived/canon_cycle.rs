use std::collections::HashSet;
use itertools::Itertools;
use yui_core::{Ring, RingOps};
use yui_link::{Link, LinkComp, State};

use crate::{KhAlgGen, KhAlgStr, KhComplex, KhLabel, KhEnhState, KhChain};

impl<R> KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn canon_cycle(&self) -> KhChain<R> {
        canon_cycle(self.link(), self.structure(), true)
    }

    pub fn canon_cycles(&self) -> Vec<KhChain<R>> {
        let ori = if self.is_reduced() { 
            vec![true]
        } else { 
            vec![true, false]
        };

        ori.into_iter().map(|o| { 
            canon_cycle(self.link(), self.structure(), o)
        }).collect()
    }
}

// -- private -- //

fn canon_cycle<R>(link: &Link, str: &KhAlgStr<R>, positive: bool) -> KhChain<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    assert!(str.t().is_zero(), "Only supported for t = 0.");
    assert_eq!(link.components().len(), 1, "Only knots are supported.");

    let s = link.ori_pres_state();
    let circles = link.seifert_circles();
    let colors = color_circles(link, &circles, positive);

    make_chain(str, &s, &colors)
}

fn color_circles(l: &Link, circles: &Vec<LinkComp>, positive: bool) -> Vec<Color> { 
    let n = circles.len();
    let mut colors = vec![Color::A; n];

    if n == 0 {
        return colors
    }

    let mut queue = vec![];
    let mut remain: HashSet<_> = (0..n).collect();

    let e = l.first_edge().unwrap();
    let i = circles.iter().find_position(|c| c.edges().contains(e)).unwrap().0;

    queue.push(i);
    colors[i] = if positive { Color::A } else { Color::B };

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

fn make_chain<R>(str: &KhAlgStr<R>, s: &State, colors: &Vec<Color>) -> KhChain<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    let x = KhEnhState::new(s.clone(), KhLabel::empty());
    let mut z = KhChain::from((x, R::one()));
    for c in colors {
        z *= c.as_tensor_factor(str);
    }

    z
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

        fn init(x: KhAlgGen) -> KhEnhState { 
            let mut z = KhEnhState::init();
            z.label.push(x);
            z
        }

        match self { 
            Color::A => KhChain::from( 
                (init(X), R::one()) 
            ),
            Color::B => KhChain::from_iter([
                (init(X), R::one()), 
                (init(I), -s.h())
            ])
        }
    }
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