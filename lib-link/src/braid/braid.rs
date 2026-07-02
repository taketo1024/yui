use std::ops::{MulAssign, Mul};
use auto_impl_ops::auto_ops;
use delegate::delegate;
use derive_more::{Display, Debug};
use itertools::Itertools;

use crate::{Link, LinkBuilder, NodeType};

use super::braid_gen::{BraidGen, from_raw};

#[derive(Clone, PartialEq, Eq, Display, Debug)]
#[display("{:?}", elements)]
pub struct Braid {
    strands: usize,
    elements: Vec<BraidGen>
}

impl Braid {
    pub fn new(strands: usize, elements: Vec<BraidGen>) -> Self {
        Self { strands, elements }
    }

    pub fn id(strands: usize) -> Self {
        Self::new(
            strands,
            vec![]
        )
    }

    pub fn generator(strands: usize, index: usize) -> Self {
        let val = i8::try_from(index).expect("index must fit in i8");
        Self::new(strands, vec![BraidGen::new(val)])
    }

    pub fn strands(&self) -> usize {
        self.strands
    }

    pub fn elements(&self) -> &[BraidGen] {
        &self.elements
    }

    delegate! {
        to self.elements {
            pub fn len(&self) -> usize;
            #[call(is_empty)]
            pub fn is_id(&self) -> bool;
        }
    }

    pub fn inv(&self) -> Self {
        Self::new(
            self.strands,
            self.elements.iter().rev().map(
                |g| g.inv()
            ).collect()
        )
    }

    pub fn extend(&self, by: usize) -> Self {
        Self::new(self.strands + by, self.elements.clone())
    }

    pub fn reduced(&self) -> Self {
        let mut stack: Vec<BraidGen> = Vec::new();
        for g in self.elements.iter().copied() {
            match stack.last() {
                Some(&top) if top.inv() == g => { stack.pop(); }
                _ => stack.push(g),
            }
        }
        Self::new(self.strands, stack)
    }

    // Plat-free braid closure: one crossing per generator; each strand position is closed by wiring
    // its visits cyclically (the wrap-around edge is the closure arc, so a single visit wires to
    // itself), and an unvisited position becomes a free loop. Crossing ports are SW=0, SE=1, NE=2, NW=3.
    pub fn closure(&self) -> Link {
        let mut b = LinkBuilder::new();

        let xs: Vec<_> = self.elements.iter().map(|g| {
            let nt = if g.sign().is_positive() {
                NodeType::XR
            } else {
                NodeType::XL
            };
            b.add_crossing(nt)
        }).collect();

        // the (top, bottom) ports each crossing presents to its two strand positions:
        // (NW = 3, SW = 0) at position i, (NE = 2, SE = 1) at position i + 1.
        let strands = self.elements.iter().zip(&xs).fold(
            vec![vec![]; self.strands],
            |mut strands, (g, &x)| {
                let i = g.index() - 1;
                strands[i]    .push(((x, 3), (x, 0)));
                strands[i + 1].push(((x, 2), (x, 1)));
                strands
            }
        );

        for visits in strands {
            if visits.is_empty() {
                b.add_loop();
            } else {
                for (&(_, bot), &(top, _)) in visits.iter().circular_tuple_windows() {
                    b.connect(bot, top);
                }
            }
        }

        // canonical braid orientation: strands run downward, entering every crossing from the
        // top ports (NW = 3, NE = 2).
        b.build_with(|_, j| j >= 2).unwrap()
    }

    pub fn display(&self) -> String {
        fn row(strands: usize, g: &BraidGen) -> String {
            let index = g.index();
            let sign = g.sign();

            (0..3).map(|r| {
                (1..=strands).map(|i| {
                    if i == index {
                        match r {
                            0 => "\\ /",
                            1 => if sign.is_positive() { " / " } else { " \\ " },
                            _ => "/ \\",
                        }
                    } else if i == index + 1 {
                        " "
                    } else {
                        "| "
                    }
                }).join("")
            }).join("\n")
        }

        self.elements.iter().map(|g|
            row(self.strands, g)
        ).join("\n")
    }

    pub fn load(name: &str) -> Result<Braid, Box<dyn std::error::Error>> {
        let json = yui_core::util::data_dir::load_json("braid", name)?;
        let code: Vec<i32> = serde_json::from_str(&json)?;
        Ok(Braid::from_iter(code))
    }
}

macro_rules! impl_from_int {
    ($($t:ty),* $(,)?) => {
        $(
            impl<const N: usize> From<[$t; N]> for Braid {
                fn from(value: [$t; N]) -> Self {
                    Self::from_iter(value)
                }
            }

            impl FromIterator<$t> for Braid {
                fn from_iter<T: IntoIterator<Item = $t>>(iter: T) -> Self {
                    let elements = iter.into_iter().map(|v|
                        from_raw(i8::try_from(v).expect("BraidGen value must fit in i8"))
                    ).collect_vec();
                    let strands = elements.iter().map(|g| g.index() + 1).max().unwrap_or(0);
                    Self::new(strands, elements)
                }
            }
        )*
    };
}

impl_from_int!(i8, i16, i32, i64);

#[auto_ops]
impl MulAssign<&Braid> for Braid {
    fn mul_assign(&mut self, rhs: &Braid) {
        assert_eq!(self.strands, rhs.strands);
        self.elements.extend(rhs.elements.iter().cloned());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use yui_core::poly::LPoly;
    use crate::misc::jones_polynomial;

    type P = LPoly<'q', i32>;

    #[test]
    fn closure_orientation() {
        // σ₁⁻⁴ closes to the (2,4) torus link with parallel downward strands. jones detects the
        // component orientations, so this locks the closure's convention (verified identical to the
        // pre-LinkBuilder closure; L4a1's PD orientation is the anti-parallel variant, NOT this one).
        let l = Braid::from([-1, -1, -1, -1]).closure();
        assert_eq!(l.n_comps(), 2);
        assert!(l.is_oriented());
        assert_eq!(l.writhe(), -4);

        let q = P::mono;
        assert_eq!(jones_polynomial(&l), P::from_iter([(q(-2), 1), (q(-4), 1), (q(-6), 1), (q(-12), 1)]));
    }

    #[test]
    fn init_by_code() {
        let b = Braid::from([1, 1, -2, -1, 3]);
        assert_eq!(b.strands(), 4);
        assert_eq!(b.len(), 5);
    }

    #[test]
    fn to_string() {
        let b = Braid::from([1, 1, -2, -1, 3]);
        assert_eq!(b.to_string(), "[1, 1, -2, -1, 3]");
    }

    #[test]
    fn display() {
        let b = Braid::from([1, 1, -2, -1, 3]);
        let display = b.display();
        assert_ne!(display, "")
    }

    #[test]
    fn reduced_empty() {
        let b = Braid::from([] as [i32; 0]);
        assert_eq!(b.reduced(), b);
    }

    #[test]
    fn reduced_no_cancel() {
        let b = Braid::from([1, 2, 3]);
        assert_eq!(b.reduced(), b);
    }

    #[test]
    fn reduced_single_pair() {
        let b = Braid::from([1, -1]);
        let r = b.reduced();
        assert!(r.is_id());
        assert_eq!(r.strands(), b.strands());
    }

    #[test]
    fn reduced_cascading() {
        let b = Braid::from([1, 2, -2, -1]);
        assert!(b.reduced().is_id());
    }

    #[test]
    fn reduced_mid_sequence() {
        let b = Braid::from([1, 2, -2, 3]);
        assert_eq!(b.reduced(), Braid::new(b.strands(), vec![BraidGen::new(1), BraidGen::new(3)]));
    }

    #[test]
    fn reduced_non_inverse_pair() {
        // σ_1 σ_2^{-1} does not reduce (different indices).
        let b = Braid::from([1, -2]);
        assert_eq!(b.reduced(), b);
    }

    #[test]
    fn closure() {
        let b = Braid::test_data("3_1");
        let l = b.closure();

        assert_eq!(l.n_crossings(), 3);
        assert_eq!(l.writhe(), 3);
        assert_eq!(l.n_comps(), 1);
    }

    #[test]
    fn extend() {
        let b = Braid::from([1, 2]).extend(2);
        assert_eq!(b.strands(), 5);
        assert_eq!(b.len(), 2);
    }

    #[test]
    fn extend_zero() {
        let b0 = Braid::from([1, 2]);
        let b1 = b0.extend(0);
        assert_eq!(b1, b0);
    }

    #[test]
    fn extend_closure_adds_loops() {
        // extending a braid by k adds k free loops to the closure.
        let b = Braid::from([1, 1, 1]).extend(2);
        let l = b.closure();
        assert_eq!(l.n_crossings(), 3);
        assert_eq!(l.n_loops(), 2);
        assert_eq!(l.n_comps(), 1 + 2); // original trefoil component + 2 free loops
    }

    #[test]
    fn closure_identity_braid() {
        // closure of the identity braid on n strands is an unlink of n components.
        let b = Braid::id(3);
        let l = b.closure();

        assert_eq!(l.n_crossings(), 0);
        assert_eq!(l.n_loops(), 3);
        assert_eq!(l.n_comps(), 3);
    }

    #[test]
    fn closure_with_free_strand() {
        // 3 strands, σ_1 once: strands 1-2 form the σ_1 closure (one component),
        // strand 3 is a free loop.
        let b = Braid::new(3, vec![BraidGen::new(1)]);
        let l = b.closure();

        assert_eq!(l.n_crossings(), 1);
        assert_eq!(l.n_loops(), 1);
        assert_eq!(l.n_comps(), 2);
    }
}
