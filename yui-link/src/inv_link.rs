use std::collections::HashMap;

use itertools::Itertools;
use crate::{Edge, Link, XCode};

// Involutive link
#[derive(Debug, Clone)]
pub struct InvLink { 
    link: Link,
    base_pt: Option<Edge>,
    e_map: HashMap<Edge, Edge>,
    x_map: HashMap<usize, usize>
}

impl InvLink { 
    pub fn new<I>(link: Link, symm: I, base_pt: Option<Edge>) -> InvLink
    where I: IntoIterator<Item = (Edge, Edge)> { 
        let mut edges = link.edges();
        let mut e_map = HashMap::new();

        for (e1, e2) in symm { 
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

    pub fn from_code<I1, I2>(pd_code: I1, symm: I2, base_pt: Option<Edge>) -> Self
    where I1: IntoIterator<Item = XCode>, I2: IntoIterator<Item = (Edge, Edge)> { 
        let l = Link::from_pd_code(pd_code);
        Self::new(l, symm, base_pt)
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
            "3_1" => Ok(InvLink::from_code(
                [[1,5,2,4],[3,1,4,6],[5,3,6,2]],
                [(2,6),(3,5)],
                Some(1)
            )),
            "4_1" => Ok(InvLink::from_code(
                [[2,7,3,8],[4,2,5,1],[6,3,7,4],[8,6,1,5]],
                [(2,8),(3,7),(4,6)],
                Some(1)
            )),
            "5_1" => Ok(InvLink::from_code(
                [[1,7,2,6],[3,9,4,8],[5,1,6,10],[7,3,8,2],[9,5,10,4]],
                [(2,10),(3,9),(4,8),(5,7)],
                Some(1)
            )), 
            "5_2a" => Ok(InvLink::from_code(
                [[3,11,4,10],[5,9,6,8],[6,2,7,1],[9,5,10,4],[11,3,12,2],[12,8,1,7]],
                [(2,12),(3,11),(4,10),(5,9),(6,8)],
                Some(1)
            )), 
            "5_2b" => Ok(InvLink::from_code(
                [[1,7,2,6],[4,10,5,9],[5,3,6,2],[7,1,8,12],[10,4,11,3],[11,9,12,8]],
                [(2,12),(3,11),(4,10),(5,9),(6,8)],
                Some(1)
            )), 
            "6_1a" => Ok(InvLink::from_code(
                [[1,6,2,7],[3,11,4,10],[5,9,6,8],[7,12,8,1],[9,5,10,4],[11,3,12,2]],
                [(2,12),(3,11),(4,10),(5,9),(6,8)],
                Some(1)
            )), 
            "6_1b" => Ok(InvLink::from_code(
                [[1,7,2,6],[3,10,4,11],[5,3,6,2],[7,1,8,12],[9,4,10,5],[11,9,12,8]],
                [(2,12),(3,11),(4,10),(5,9),(6,8)],
                Some(1)
            )), 
            "6_2a" => Ok(InvLink::from_code(
                [[1,6,2,7],[3,11,4,10],[5,9,6,8],[7,12,8,1],[9,3,10,2],[11,5,12,4]],
                [(2,12),(3,11),(4,10),(5,9),(6,8)],
                Some(1)
            )), 
            "6_2b" => Ok(InvLink::from_code(
                [[1,9,2,8],[4,11,5,12],[7,1,8,14],[9,3,10,2],[10,5,11,6],[12,3,13,4],[13,7,14,6]],
                [(2,14),(3,13),(4,12),(5,11),(6,10),(7,9)],
                Some(1)
            )), 
            "6_3" => Ok(InvLink::from_code(
                [[3,13,4,12],[6,9,7,10],[8,1,9,2],[10,5,11,6],[11,3,12,2],[13,5,14,4],[14,7,1,8]],
                [(2,14),(3,13),(4,12),(5,11),(6,10),(7,9)],
                Some(1)
            )), 
            "7_1" => Ok(InvLink::from_code(
                [[1,9,2,8],[3,11,4,10],[5,13,6,12],[7,1,8,14],[9,3,10,2],[11,5,12,4],[13,7,14,6]],
                [(2,14),(3,13),(4,12),(5,11),(6,10),(7,9)],
                Some(1)
            )), 
            "7_2a" => Ok(InvLink::from_code(
                [[3,15,4,14],[5,13,6,12],[7,11,8,10],[8,2,9,1],[11,7,12,6],[13,5,14,4],[15,3,16,2],[16,10,1,9]],
                [(2,16),(3,15),(4,14),(5,13),(6,12),(7,11),(8,10)],
                Some(1)
            )), 
            "7_2b" => Ok(InvLink::from_code(
                [[1,9,2,8],[3,7,4,6],[4,14,5,13],[7,3,8,2],[9,1,10,16],[11,15,12,14],[12,6,13,5],[15,11,16,10]],
                [(2,16),(3,15),(4,14),(5,13),(6,12),(7,11),(8,10)],
                Some(1)
            )), 
            "7_3a" => Ok(InvLink::from_code(
                [[1,9,2,8],[3,13,4,12],[5,11,6,10],[7,1,8,14],[9,3,10,2],[11,5,12,4],[13,7,14,6]],
                [(2,14),(3,13),(4,12),(5,11),(6,10),(7,9)],
                Some(1)
            )), 
            "7_3b" => Ok(InvLink::from_code(
                [[3,13,4,12],[5,15,6,14],[8,2,9,1],[10,7,11,8],[11,3,12,2],[13,5,14,4],[15,7,16,6],[16,10,1,9]],
                [(2,16),(3,15),(4,14),(5,13),(6,12),(7,11),(8,10)],
                Some(1)
            )), 
            "7_4a" => Ok(InvLink::from_code(
                [[2,8,3,7],[3,15,4,14],[5,13,6,12],[8,2,9,1],[10,16,11,15],[11,7,12,6],[13,5,14,4],[16,10,1,9]],
                [(2,16),(3,15),(4,14),(5,13),(6,12),(7,11),(8,10)],
                Some(1)
            )), 
            "7_4b" => Ok(InvLink::from_code(
                [[2,10,3,9],[4,12,5,11],[6,14,7,13],[8,4,9,3],[10,2,11,1],[12,8,13,7],[14,6,1,5]],
                [(2,14),(3,13),(4,12),(5,11),(6,10),(7,9)],
                Some(1)
            )), 
            "7_5a" => Ok(InvLink::from_code(
                [[1,9,2,8],[3,13,4,12],[5,11,6,10],[7,1,8,14],[9,7,10,6],[11,3,12,2],[13,5,14,4]],
                [(2,14),(3,13),(4,12),(5,11),(6,10),(7,9)],
                Some(1)
            )), 
            "7_5b" => Ok(InvLink::from_code(
                [[1,11,2,10],[4,8,5,7],[5,15,6,14],[9,1,10,18],[11,3,12,2],[12,16,13,15],[13,7,14,6],[16,3,17,4],[17,9,18,8]],
                [(2,18),(3,17),(4,16),(5,15),(6,14),(7,13),(8,12),(9,11)],
                Some(1)
            )), 
            "7_6a" => Ok(InvLink::from_code(
                [[2,13,3,14],[4,11,5,12],[6,4,7,3],[8,1,9,2],[10,5,11,6],[12,10,13,9],[14,7,1,8]],
                [(2,14),(3,13),(4,12),(5,11),(6,10),(7,9)],
                Some(1)
            )), 
            "7_6b" => Ok(InvLink::from_code(
                [[1,8,2,9],[3,15,4,14],[6,11,7,12],[7,4,8,5],[9,16,10,1],[12,5,13,6],[13,10,14,11],[15,3,16,2]],
                [(2,16),(3,15),(4,14),(5,13),(6,12),(7,11),(8,10)],
                Some(1)
            )), 
            "7_7a" => Ok(InvLink::from_code(
                [[1,8,2,9],[4,13,5,14],[7,11,8,10],[9,16,10,1],[11,3,12,2],[12,5,13,6],[14,3,15,4],[15,7,16,6]],
                [(2,16),(3,15),(4,14),(5,13),(6,12),(7,11),(8,10)],
                Some(1)
            )), 
            "7_7b" => Ok(InvLink::from_code(
                [[1,10,2,11],[3,13,4,12],[5,14,6,1],[7,5,8,4],[9,2,10,3],[11,9,12,8],[13,6,14,7]],
                [(2,14),(3,13),(4,12),(5,11),(6,10),(7,9)],
                Some(1)
            )), 
            _ => todo!()
        }
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