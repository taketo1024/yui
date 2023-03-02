use crate::{Edge, Link};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Component { 
    edges:Vec<Edge>,
    closed:bool
}

impl Component { 
    pub fn new(edges: Vec<Edge>, closed: bool) -> Self { 
        Self { edges, closed }
    }

    pub fn edges(&self) -> &Vec<Edge> { 
        &self.edges
    }

    pub fn is_adj(&self, other: &Component, link: &Link) -> bool { 
        // 1) find crossings `x` that touche `self`. 
        // 2) check if `x` also touches `other`.
        
        for x in link.data().iter() { 
            if !x.edges().iter().any(|e| 
                self.edges.contains(e)
            ) { 
                continue
            }

            let Some(e) = x.edges().iter().find(|e| 
                !self.edges.contains(e)
            ) else { 
                continue
            };

            if other.edges.contains(e) { 
                return true
            }
        }

        false
    }

    pub fn is_closed(&self) -> bool { 
        self.closed
    }
}