use std::collections::HashMap;
use std::fmt;
use std::ops::{MulAssign, Mul};
use auto_impl_ops::auto_ops;
use delegate::delegate;
use itertools::Itertools;
use num_traits::Zero;
use yui::{GetSign, Sign};
use yui::format::{subscript, superscript};

use crate::{Link, XCode};

#[derive(Clone, Copy, PartialEq, Eq)]
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
        self.0.abs() as usize
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

impl fmt::Display for Generator {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.sign().is_positive() {
            write!(f, "σ{}", subscript(self.index()))
        } else { 
            write!(f, "σ{}{}", subscript(self.index()), superscript(-1))
        }
    }
}

impl fmt::Debug for Generator {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(&self, f)
    }
}

#[derive(Clone, PartialEq, Eq)]
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

    pub fn gen(strands: usize, index: usize) -> Self {
        Self::new(
            strands,
            vec![(Generator::new(index, Sign::Pos))]
        )
    }

    pub fn all_gens(strands: usize) -> Vec<Self> {
        (1..strands).map(|index| Self::gen(strands, index)).collect()
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
                |gen| gen.inv()
            ).collect()
        )
    }

    pub fn reduce(&mut self) {
        // TODO
    }

    pub fn closure(&self) -> Link {
        let mut count = self.strands;
        let mut bottom_edges: Vec<usize> = (0..self.strands).collect();
        let mut pd_code: Vec<XCode> = Vec::new();

        for s in &self.elements {
            /*        +       -
             *  ↓   a   b   a   b
             *       \ /     \ /
             *  ↓     /       \
             *       / \     / \
             *  ↓   c   d   c   d
             */
            let i = s.index() - 1;
            let (a, b) = (bottom_edges[i], bottom_edges[i + 1]);
            let (c, d) = (count, count + 1);

            if s.sign().is_positive() {
                pd_code.push([a, c, d, b]);
            } else {
                pd_code.push([b, a, c, d]);
            }

            bottom_edges[i] = c;
            bottom_edges[i + 1] = d;
            count += 2;
        }

        assert!(
            bottom_edges.iter().enumerate().all(|(i, &j)| i != j),
            "braid closure contains free loop."
        );

        let conn: HashMap<_, _> = Iterator::zip(
            bottom_edges.into_iter(),
            0..self.strands
        ).collect();

        let pd_code = pd_code
            .into_iter()
            .map(|x| x.map(|a| *conn.get(&a).unwrap_or(&a)))
            .collect_vec();

        Link::from_pd_code(pd_code)
    }

    pub fn display(&self) -> String { 
        fn row(strands: usize, gen: &Generator) -> String {
            let index = gen.index();
            let sign = gen.sign();

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

        self.elements.iter().map(|gen| 
            row(self.strands, gen)
        ).join("\n")
    }

    pub fn load(name_or_path: &str) -> Result<Braid, Box<dyn std::error::Error>> {
        const RESOURCE_DIR: &str = "resources/braid/";
        
        use regex::Regex;
        let r = Regex::new(r"([1-9]|10)_[0-9]+").unwrap();

        if r.is_match(name_or_path) { 
            let dir = std::env!("CARGO_MANIFEST_DIR");
            let path = format!("{dir}/{RESOURCE_DIR}{name_or_path}.json");
            Self::_load(&path)
        } else { 
            Self::_load(name_or_path)
        }
    }

    fn _load(path: &str) -> Result<Braid, Box<dyn std::error::Error>> {
        let json = std::fs::read_to_string(path)?;
        let code: Vec<i32> = serde_json::from_str(&json)?;
        let braid = Braid::from_iter(code);
        Ok(braid)
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

impl fmt::Display for Braid {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for gen in self.elements.iter() { // should reverse?
            write!(f, "{gen}")?
        }
        Ok(())
    }
}

impl fmt::Debug for Braid {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(&self, f)
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
        println!("{}", b.display());
    }

    #[test]
    fn load() {
        let name = "3_1";
        let b = Braid::load(name);

        assert!(b.is_ok());

        let b = b.unwrap();
        assert_eq!(b.strands(), 2);
    }

    #[test]
    fn closure() {
        let b = Braid::from([-1,-1,-2,1,3,2,2,-4,-3,2,-3,-4]); // 9_41
        let l = b.closure();
        assert_eq!(l.crossing_num(), 12);
    }
}