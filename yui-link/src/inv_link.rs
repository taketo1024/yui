use std::collections::HashMap;

use itertools::Itertools;
use crate::{Edge, Link};

// Involutive link
#[derive(Debug, Clone)]
pub struct InvLink { 
    link: Link,
    base_pt: Option<Edge>,
    e_map: HashMap<Edge, Edge>,
    x_map: HashMap<usize, usize>
}

impl InvLink { 
    pub fn new<I>(link: Link, e_pairs: I, base_pt: Option<Edge>) -> InvLink
    where I: IntoIterator<Item = (Edge, Edge)> { 
        let mut edges = link.edges();
        let mut e_map = HashMap::new();

        for (e1, e2) in e_pairs { 
            assert!(edges.contains(&e1));
            assert!(edges.contains(&e2));

            e_map.insert(e1, e2);
            e_map.insert(e2, e1);
            edges.remove(&e1);
            edges.remove(&e2);
        }

        for e in edges { 
            e_map.insert(e, e);
        }

        let mut x_map = HashMap::new();

        // TODO? check resolution
        for (i, x) in link.data().iter().enumerate() { 
            let edges = x.edges().map(|e| e_map.get(&e).unwrap());
            let find = link.data().iter().find_position(|y|
                edges.iter().all(|e| y.edges().contains(e))
            );

            assert!(find.is_some(), "no match for x: {x} -> {edges:?}");

            let j = find.unwrap().0;
            x_map.insert(i, j);
            x_map.insert(j, i);
        }

        Self { link, base_pt, e_map, x_map }
    }

    pub fn link(&self) -> &Link { 
        &self.link
    }

    pub fn base_pt(&self) -> Option<Edge> {
        self.base_pt
    }

    pub fn inv_e(&self, e: Edge) -> Edge { 
        self.e_map.get(&e).cloned().unwrap()
    }

    pub fn inv_x(&self, i: usize) -> usize { 
        self.x_map.get(&i).cloned().unwrap()
    }

    pub fn mirror(&self) -> Self { 
        Self { 
            link: self.link.mirror(),
            base_pt: self.base_pt,
            e_map: self.e_map.clone(),
            x_map: self.x_map.clone(), 
        }
    }
}

impl InvLink {
    pub fn load(name: &str) -> Result<InvLink, Box<dyn std::error::Error>> {
        match name {
            _ => todo!()
        }
    }
}

impl InvLink { 
    #[allow(non_snake_case)]
    pub fn knotJ() -> Self { 
        InvLink::new(
            Link::from_pd_code([
                [0,26,1,25],[18,1,19,2],[2,12,3,11],[3,30,4,31],
                [29,4,30,5],[12,6,13,5],[7,26,8,27],[8,0,9,33],
                [9,17,10,16],[23,10,24,11],[13,20,14,21],[27,15,28,14],
                [32,15,33,16],[17,25,18,24],[19,7,20,6],[28,22,29,21],
                [22,32,23,31]
            ]), 
            [(1,33),(2,32),(3,31),(4,30),(5,29),(6,28),(7,27),(8,26),(9,25),(10,24),(11,23),(12,22),(13,21),(14,20),(15,19),(16,18)],
            Some(0)
        )
    }
}

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn inv_e() { 
        let l = Link::from_pd_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]);
        let l = InvLink::new(l, [(1,5), (2,4)], None);

        assert_eq!(l.inv_e(1), 5);
        assert_eq!(l.inv_e(2), 4);
        assert_eq!(l.inv_e(3), 3);
        assert_eq!(l.inv_e(4), 2);
        assert_eq!(l.inv_e(5), 1);
        assert_eq!(l.inv_e(6), 6);
    }
    
    #[test]
    fn inv_x() { 
        let l = Link::from_pd_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]);
        let l = InvLink::new(l, [(1,5), (2,4)], None);

        assert_eq!(l.inv_x(0), 0);
        assert_eq!(l.inv_x(1), 2);
        assert_eq!(l.inv_x(2), 1);
    }
}