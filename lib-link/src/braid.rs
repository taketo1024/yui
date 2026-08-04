use std::collections::HashMap;
use std::ops::{MulAssign, Mul};
use auto_impl_ops::auto_ops;
use delegate::delegate;
use derive_more::{Display, Debug};
use itertools::Itertools;
use num_traits::Zero;
use yui_core::{GetSign, Sign};

use crate::{Link, Node};

#[derive(Clone, Copy, PartialEq, Eq, Display, Debug)]
#[display("{}", _0)]
#[  debug("{}", _0)]
pub struct Generator(i32);

impl Generator { 
    pub fn new(index: usize, sign: Sign) -> Self { 
        assert!(!index.is_zero());
        if sign.is_positive() { 
            Self(index as i32 )
        } else { 
            Self(-(index as i32))
        }
    }

    pub fn index(&self) -> usize {
        self.0.unsigned_abs() as usize
    }
    
    pub fn sign(&self) -> Sign { 
        self.0.sign()
    }

    pub fn inv(&self) -> Self { 
        Self(-self.0)
    }
}

impl From<i32> for Generator {
    fn from(value: i32) -> Self {
        assert!(!value.is_zero());
        Self(value)
    }
}

#[derive(Clone, PartialEq, Eq, Display, Debug)]
#[display("{:?}", elements)]
pub struct Braid {
    strands: usize,
    elements: Vec<Generator>
}

impl Braid {
    pub fn new(strands: usize, elements: Vec<Generator>) -> Self {
        Self { strands, elements }
    }

    pub fn id(strands: usize) -> Self {
        Self::new(
            strands,
            vec![]
        )
    }

    pub fn generator(strands: usize, index: usize) -> Self {
        Self::new(
            strands,
            vec![(Generator::new(index, Sign::Pos))]
        )
    }

    pub fn strands(&self) -> usize { 
        self.strands
    }

    pub fn elements(&self) -> &[Generator] { 
        &self.elements
    }

    delegate! { 
        to self.elements { 
            pub fn len(&self) -> usize;
            #[call(is_empty)] 
            pub fn is_triv(&self) -> bool;
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

    pub fn reduce(&mut self) {
        // TODO
    }

    pub fn closure(&self) -> Link {
        use crate::{NodeType, NodeOri};

        let mut count = self.strands;
        let mut front_edges: Vec<usize> = (0..self.strands).collect();
        let mut nodes: Vec<Node> = Vec::new();

        for s in &self.elements {
            /*        +       -
             *  ↓   a   b   a   b
             *       \ /     \ /
             *  ↓     /       \
             *       / \     / \
             *  ↓   c   d   c   d
             */
            let i = s.index() - 1;
            let (a, b) = (front_edges[i], front_edges[i + 1]);
            let (c, d) = (count, count + 1);

            let nt = if s.sign().is_positive() { NodeType::XR } else { NodeType::XL };
            let n = Node::new(nt, NodeOri::Up, [b,a,c,d]);
            nodes.push(n);

            front_edges[i] = c;
            front_edges[i + 1] = d;
            count += 2;
        }

        assert!(
            front_edges.iter().enumerate().all(|(i, &j)| i != j),
            "braid closure contains free loop."
        );

        let conn: HashMap<_, _> = Iterator::zip(
            front_edges.into_iter(),
            0..self.strands
        ).collect();

        nodes.iter_mut().for_each(|n|
            *n = n.convert_edges(|e| conn.get(&e).cloned().unwrap_or(e))
        );

        Link::from_nodes(nodes)
    }

    pub fn display(&self) -> String { 
        fn row(strands: usize, g: &Generator) -> String {
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

impl<const N: usize> From<[i32; N]> for Braid {
    fn from(value: [i32; N]) -> Self {
        Self::from_iter(value)
    }
}

impl FromIterator<i32> for Braid {
    fn from_iter<T: IntoIterator<Item = i32>>(iter: T) -> Self {
        let elements = iter.into_iter().map(Generator).collect_vec();
        let strands = elements.iter().map(|g| g.index() + 1).max().unwrap_or(0);
        Self::new(strands, elements)
    }
}

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
    fn closure() {
        let b = Braid::test_data("3_1");
        let l = b.closure();

        assert_eq!(l.n_crossings(), 3);
        assert_eq!(l.writhe(), 3);
        assert_eq!(l.n_comps(), 1);
    }
}