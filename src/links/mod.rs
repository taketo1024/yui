use CrossingType::{Xp, Xn, V, H};
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

#[derive(Debug)]
pub struct Link { 

}

impl Link {
    pub fn new() -> Link { 
        Link{}
    }
}

#[cfg(test)]
mod tests;