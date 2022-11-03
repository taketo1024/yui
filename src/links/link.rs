
#[derive(Debug)]
struct Edge { 
    id: u8
}

#[derive(Debug)]
enum Resolution { 
    Res0, Res1
}

#[derive(Debug)]
struct State { 
    res: Vec<Resolution>
}

#[derive(Debug, Clone, Copy)]
enum CrossingType { 
    Xp, Xn, V, H 
}

impl CrossingType { 
    fn mirror(self) -> CrossingType {
        use CrossingType::{Xp, Xn};
        match self { 
            Xp => Xn,
            Xn => Xp, 
            other => other
        }
    }
}

struct Crossing { 
    ctype: CrossingType,
    edges: (Edge, Edge, Edge, Edge)
}

impl Crossing {
    fn is_resolved(&self) -> bool { 
        use CrossingType::{V, H};
        match self.ctype { 
            V | H => true,
            _ => false
        }
    }

    fn resolve(&mut self, r: Resolution) {
        use CrossingType::{Xp, Xn, V, H};
        use Resolution::{Res0, Res1};
        
        match (self.ctype, r) {
            (Xp, Res0) | (Xn, Res1) => self.ctype = V,
            (Xp, Res1) | (Xn, Res0) => self.ctype = H,
            _ => panic!()
        }
    }

    fn mirror(&mut self) { 
        self.ctype = self.ctype.mirror();
    }

    fn pass(index:u8) -> u8 { 
        fn p(j: u8) -> u8 { 
            0
        }
        p(index - 1) + index
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