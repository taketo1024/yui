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

    fn pass(&self, index:usize) -> usize { 
        debug_assert!((0..4).contains(&index));

        match self.ctype {
            Xp | Xn => (index + 2) % 4,
            V => 3 - index,
            H => (5 - index) % 4
        }
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

    pub fn writhe(&self) -> i8 { 
        self.crossing_signs().iter().sum()
    }

    fn next(&self, c_index:usize, e_index:usize) -> Option<(usize, usize)> {
        let n = self.data.len();
        debug_assert!((0..n).contains(&c_index));
        debug_assert!((0..4).contains(&e_index));

        let e = &self.data[c_index].edges[e_index];

        for (i, c) in self.data.iter().enumerate() { 
            for (j, f) in c.edges.iter().enumerate() { 
                if e == f && (c_index != i || (c_index == i && e_index != j)) { 
                    return Some((i, j))
                }
            }
        }
        None
    }

    fn traverse_edges<F>(&self, c_index:usize, e_index:usize, mut f:F) where
        F: FnMut(usize, usize)
    {
        let n = self.data.len();
        debug_assert!((0..n).contains(&c_index));
        debug_assert!((0..4).contains(&e_index));

        let mut i = c_index;
        let mut j = e_index;

        loop {
            f(i, j);

            let c = &self.data[i];
            let k = c.pass(j);
            
            if let Some(next) = self.next(i, k) { 
                if next != (c_index, e_index) {
                    (i, j) = next;
                    continue
                } else {
                    f(c_index, e_index); // final call
                }
            }
            break
        }
    }

    fn crossing_signs(&self) -> Vec<i8> {
        let n = self.data.len();
        let mut signs = vec![0; n];

        let mut traverse = |e0: usize| {
            for i0 in 0..n {
                if signs[i0] != 0 { 
                    continue 
                }
                self.traverse_edges(i0, e0, |i, j| { 
                    let c = &self.data[i];
                    let sign = match (c.ctype, j) { 
                        (Xp, 1) | (Xn, 3) =>  1,
                        (Xp, 3) | (Xn, 1) => -1,
                        _ => 0
                    };
                    if sign != 0 {
                        signs[i] = sign;
                    }
                })
            }
        };

        traverse(0);
        traverse(1); // in case 

        signs
    }
}

#[cfg(test)]
mod tests;