use std::collections::HashMap;

use itertools::Itertools;
use num_integer::Integer;
use crate::{Crossing, Edge, Link, XCode};

// Involutive link
#[derive(Debug, Clone)]
pub struct InvLink { 
    link: Link,
    base_pt: Option<Edge>,
    e_map: HashMap<Edge, Edge>,
    x_map: HashMap<Crossing, Crossing>
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
        for x in link.data().iter() { 
            let edges = x.edges().map(|e| e_map.get(&e).unwrap());
            let find = link.data().iter().find_position(|y|
                edges.iter().all(|e| y.edges().contains(e))
            );

            assert!(find.is_some(), "no match for x: {x} -> {edges:?}");

            let j = find.unwrap().0;
            let y = &link.data()[j];

            x_map.insert(x.clone(), y.clone());

            if x != y { 
                x_map.insert(y.clone(), x.clone());
            }
        }

        assert_eq!(x_map.len(), link.data().len());

        if let Some(p) = base_pt { 
            assert_eq!(p, e_map[&p], "base-pt must be on-axis.");
        }

        Self { link, base_pt, e_map, x_map }
    }

    pub fn from_code<I1, I2>(pd_code: I1, symm: I2, base_pt: Option<Edge>) -> Self
    where I1: IntoIterator<Item = XCode>, I2: IntoIterator<Item = (Edge, Edge)> { 
        let l = Link::from_pd_code(pd_code);
        Self::new(l, symm, base_pt)
    }

    pub fn sinv_knot_from_code<I1>(pd_code: I1) -> Self
    where I1: IntoIterator<Item = XCode> { 
        let code = pd_code.into_iter().collect_vec();
        let l = Link::from_pd_code(code);
        let n = l.edges().len();

        assert!(n.is_even(), "number of edges must be even.");
        assert_eq!(l.edges().iter().min(), Some(&1), "edge must start from index 1.");
        assert_eq!(l.edges().iter().max(), Some(&n), "edges must have sequential indexing.");

        let symm = (2..=n/2).map(|i| (i, n - i + 2));
        
        Self::new(l, symm, Some(1))
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

    pub fn inv_x(&self, x: &Crossing) -> &Crossing { 
        self.x_map.get(x).unwrap()
    }

    pub fn mirror(&self) -> Self { 
        Self { 
            link: self.link.mirror(),
            base_pt: self.base_pt,
            e_map: self.e_map.clone(),
            x_map: self.x_map.iter().map(|(x, y)| (x.mirror(), y.mirror())).collect(), 
        }
    }
}

impl InvLink {
    pub fn load(name: &str) -> Result<InvLink, Box<dyn std::error::Error>> {
        match name {
            "3_1" => Ok(InvLink::sinv_knot_from_code(
                [[1,5,2,4],[3,1,4,6],[5,3,6,2]]
            )),
            "4_1" => Ok(InvLink::sinv_knot_from_code(
                [[2,7,3,8],[4,2,5,1],[6,3,7,4],[8,6,1,5]],
            )),
            "5_1" => Ok(InvLink::sinv_knot_from_code(
                [[1,7,2,6],[3,9,4,8],[5,1,6,10],[7,3,8,2],[9,5,10,4]],
            )), 
            "5_2a" => Ok(InvLink::sinv_knot_from_code(
                [[3,11,4,10],[5,9,6,8],[6,2,7,1],[9,5,10,4],[11,3,12,2],[12,8,1,7]],
            )), 
            "5_2b" => Ok(InvLink::sinv_knot_from_code(
                [[1,7,2,6],[4,10,5,9],[5,3,6,2],[7,1,8,12],[10,4,11,3],[11,9,12,8]],
            )), 
            "6_1a" => Ok(InvLink::sinv_knot_from_code(
                [[1,6,2,7],[3,11,4,10],[5,9,6,8],[7,12,8,1],[9,5,10,4],[11,3,12,2]],
            )), 
            "6_1b" => Ok(InvLink::sinv_knot_from_code(
                [[1,7,2,6],[3,10,4,11],[5,3,6,2],[7,1,8,12],[9,4,10,5],[11,9,12,8]],
            )), 
            "6_2a" => Ok(InvLink::sinv_knot_from_code(
                [[1,6,2,7],[3,11,4,10],[5,9,6,8],[7,12,8,1],[9,3,10,2],[11,5,12,4]],
            )), 
            "6_2b" => Ok(InvLink::sinv_knot_from_code(
                [[1,9,2,8],[4,11,5,12],[7,1,8,14],[9,3,10,2],[10,5,11,6],[12,3,13,4],[13,7,14,6]],
            )), 
            "6_3" => Ok(InvLink::sinv_knot_from_code(
                [[3,13,4,12],[6,9,7,10],[8,1,9,2],[10,5,11,6],[11,3,12,2],[13,5,14,4],[14,7,1,8]],
            )), 
            "7_1" => Ok(InvLink::sinv_knot_from_code(
                [[1,9,2,8],[3,11,4,10],[5,13,6,12],[7,1,8,14],[9,3,10,2],[11,5,12,4],[13,7,14,6]],
            )), 
            "7_2a" => Ok(InvLink::sinv_knot_from_code(
                [[3,15,4,14],[5,13,6,12],[7,11,8,10],[8,2,9,1],[11,7,12,6],[13,5,14,4],[15,3,16,2],[16,10,1,9]],
            )), 
            "7_2b" => Ok(InvLink::sinv_knot_from_code(
                [[1,9,2,8],[3,7,4,6],[4,14,5,13],[7,3,8,2],[9,1,10,16],[11,15,12,14],[12,6,13,5],[15,11,16,10]],
            )), 
            "7_3a" => Ok(InvLink::sinv_knot_from_code(
                [[1,9,2,8],[3,13,4,12],[5,11,6,10],[7,1,8,14],[9,3,10,2],[11,5,12,4],[13,7,14,6]],
            )), 
            "7_3b" => Ok(InvLink::sinv_knot_from_code(
                [[3,13,4,12],[5,15,6,14],[8,2,9,1],[10,7,11,8],[11,3,12,2],[13,5,14,4],[15,7,16,6],[16,10,1,9]],
            )), 
            "7_4a" => Ok(InvLink::sinv_knot_from_code(
                [[2,8,3,7],[3,15,4,14],[5,13,6,12],[8,2,9,1],[10,16,11,15],[11,7,12,6],[13,5,14,4],[16,10,1,9]],
            )), 
            "7_4b" => Ok(InvLink::sinv_knot_from_code(
                [[2,10,3,9],[4,12,5,11],[6,14,7,13],[8,4,9,3],[10,2,11,1],[12,8,13,7],[14,6,1,5]],
            )), 
            "7_5a" => Ok(InvLink::sinv_knot_from_code(
                [[1,9,2,8],[3,13,4,12],[5,11,6,10],[7,1,8,14],[9,7,10,6],[11,3,12,2],[13,5,14,4]],
            )), 
            "7_5b" => Ok(InvLink::sinv_knot_from_code(
                [[1,11,2,10],[4,8,5,7],[5,15,6,14],[9,1,10,18],[11,3,12,2],[12,16,13,15],[13,7,14,6],[16,3,17,4],[17,9,18,8]],
            )), 
            "7_6a" => Ok(InvLink::sinv_knot_from_code(
                [[2,13,3,14],[4,11,5,12],[6,4,7,3],[8,1,9,2],[10,5,11,6],[12,10,13,9],[14,7,1,8]],
            )), 
            "7_6b" => Ok(InvLink::sinv_knot_from_code(
                [[1,8,2,9],[3,15,4,14],[6,11,7,12],[7,4,8,5],[9,16,10,1],[12,5,13,6],[13,10,14,11],[15,3,16,2]],
            )), 
            "7_7a" => Ok(InvLink::sinv_knot_from_code(
                [[1,8,2,9],[4,13,5,14],[7,11,8,10],[9,16,10,1],[11,3,12,2],[12,5,13,6],[14,3,15,4],[15,7,16,6]],
            )), 
            "7_7b" => Ok(InvLink::sinv_knot_from_code(
                [[1,10,2,11],[3,13,4,12],[5,14,6,1],[7,5,8,4],[9,2,10,3],[11,9,12,8],[13,6,14,7]],
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
        let data = l.link.data();

        assert_eq!(l.inv_x(&data[0]), &data[0]);
        assert_eq!(l.inv_x(&data[1]), &data[2]);
        assert_eq!(l.inv_x(&data[2]), &data[1]);
    }

    #[test]
    fn from_sinv() { 
        let l = InvLink::sinv_knot_from_code([[1,5,2,4],[3,1,4,6],[5,3,6,2]]);

        assert_eq!(l.inv_e(1), 1);
        assert_eq!(l.inv_e(2), 6);
        assert_eq!(l.inv_e(3), 5);
        assert_eq!(l.inv_e(4), 4);
    }

    #[test]
    fn load_3_1() { 
        let l = InvLink::load("3_1").unwrap();
        assert_eq!(l.link().crossing_num(), 3);
    }

    #[test]
    fn load_4_1() { 
        let l = InvLink::load("4_1").unwrap();
        assert_eq!(l.link().crossing_num(), 4);
    }
}