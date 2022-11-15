use Resolution::{Res0, Res1};

type Edge = u8;

#[derive(Debug, PartialEq)]
enum Resolution { 
    Res0, Res1
}

#[derive(Debug, PartialEq)]
struct State { 
    res: Vec<Resolution>
}

use CrossingType::{Xp, Xn, V, H};

#[derive(Debug, Clone, Copy, PartialEq)]
enum CrossingType { 
    Xp, Xn, V, H 
}

impl CrossingType { 
    fn mirror(self) -> CrossingType {
        match self { 
            Xp => Xn,
            Xn => Xp, 
            other => other
        }
    }
}

#[derive(Debug, PartialEq)]
struct Crossing { 
    ctype: CrossingType,
    edges: [Edge; 4]
}

impl Crossing {
    fn is_resolved(&self) -> bool { 
        match self.ctype { 
            V | H => true,
            _ => false
        }
    }

    fn resolve(&mut self, r: Resolution) {
        match (self.ctype, r) {
            (Xp, Res0) | (Xn, Res1) => self.ctype = V,
            (Xp, Res1) | (Xn, Res0) => self.ctype = H,
            _ => panic!()
        }
    }

    fn mirror(&mut self) { 
        self.ctype = self.ctype.mirror();
    }

    fn pass(&self, index:usize) -> Edge { 
        debug_assert!((0..4).contains(&index));

        fn next(ctype: CrossingType, j: usize) -> usize { 
            match ctype {
                Xp | Xn => (j + 2) % 4,
                V => 3 - j,
                H => (5 - j) % 4
            }
        }

        let ctype = self.ctype;
        self.edges[ next(ctype, index) ]
    }
}

// Planer Diagram code, represented by crossings:
//
//     3   2
//      \ /
//       \      = (0, 1, 2, 3)
//      / \
//     0   1
//
// The lower edge has direction 0 -> 2.
// The crossing is +1 if the upper goes 3 -> 1.
// see: http://katlas.math.toronto.edu/wiki/Planar_Diagrams

#[derive(Debug)]
pub struct Link { 
    data: Vec<Crossing>
}

impl Link {
    pub fn empty() -> Link {
        Link { data: vec![] }
    }

    pub fn new(pd_code: &[[Edge; 4]]) -> Link { 
        let data = pd_code
            .iter()
            .map( |x| 
                Crossing { 
                    ctype: Xn,
                    edges: *x
                }
            )
            .collect();

        Link{ data }
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    pub fn crossing_num(&self) -> usize { 
        self.data.iter()
            .filter(|x| !x.is_resolved())
            .count()
    }
}

#[cfg(test)]
mod tests;