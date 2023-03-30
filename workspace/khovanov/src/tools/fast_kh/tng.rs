use std::collections::HashSet;
use std::fmt::Display;
use itertools::Itertools;
use yui_link::{Component, Edge, Crossing, Resolution};

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Tng {
    comps: Vec<Component> // arc or circle
}

impl Tng { 
    pub fn new(comps: Vec<Component>) -> Self { 
        debug_assert!(comps.iter().all(|c| !c.is_empty()));
        Self { comps }
    }

    pub fn from_x(x: &Crossing, r: Resolution) -> Self { 
        let (r0, r1) = x.res_arcs(r);
        Self::new(vec![r0, r1])
    }

    pub fn empty() -> Self { 
        Self::new(vec![])
    }

    pub fn is_empty(&self) -> bool { 
        self.comps.is_empty()
    }

    pub fn is_closed(&self) -> bool { 
        self.comps.iter().all(|a| a.is_empty() || a.is_circle())
    }

    pub fn ncomps(&self) -> usize { 
        self.comps.len()
    }

    pub fn comps(&self) -> &Vec<Component> { 
        &self.comps
    }

    pub fn comp(&self, i: usize) -> &Component { 
        &self.comps[i]
    }

    pub fn sub<I>(&self, iter: I) -> Self 
    where I: IntoIterator<Item = usize> { 
        Self::new(
            iter.into_iter().map(|i| self.comp(i).clone()).collect()
        )
    }

    pub fn endpts(&self) -> HashSet<Edge> { 
        self.comps.iter().flat_map(|c| 
            c.ends().map(|(e0, e1)| vec![e0, e1]).unwrap_or(vec![])
        ).collect()
    }

    pub fn contains(&self, c: &Component) -> bool { 
        self.comps.contains(c)
    }

    pub fn index_of(&self, c: &Component) -> Option<usize> {
        self.comps.iter().position(|c1| c1 == c)
    }

    pub fn remove_at(&mut self, i: usize) -> Component {
        self.comps.remove(i)
    }

    pub fn connect(&mut self, other: Self) {
        for c in other.comps.into_iter() {
            if c.is_circle() { 
                self.comps.push(c);
            } else { 
                self.append_arc(c);
            }
        }
    }

    pub fn append_arc(&mut self, arc: Component) { 
        assert!(arc.is_arc());

        let n = self.ncomps();

        // If one end of `arc` is connectable:
        if let Some(i) = self.find_connectable(&arc, n) { 
            self.comps[i].connect(arc);

            // If the other end is also connectable to a different component:
            if let Some(j) = self.find_connectable(&self.comps[i], i) { 
                let arc_j = self.comps.remove(j);
                self.comps[i].connect(arc_j);
            }
        } else { 
            self.comps.push(arc);
        }
    }

    fn find_connectable(&self, arc: &Component, j: usize) -> Option<usize> {
        let n = self.comps.len();

        (0..n).find(|&i| 
            i != j && self.comps[i].is_connectable(arc)
        )
    }

    pub fn find_loop(&self) -> Option<usize> {
        self.comps.iter().enumerate().find(|(_, c)| 
            c.is_circle()
        ).map(|(i, _)|
            i
        )
    }

    pub fn euler_num(&self) -> isize { 
        // NOTE: χ(arc) = 1, χ(circle) = 0. 
        self.comps.iter().filter(|c| 
            c.is_arc()
        ).count() as isize
    }
}

impl Display for Tng {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_empty() { 
            write!(f, "∅")
        } else { 
            let arcs = self.comps.iter().map(|a| 
                a.to_string()
            ).join(", ");
            write!(f, "{{{}}}", arcs)
        }
    }
}

impl From<Component> for Tng {
    fn from(c: Component) -> Self {
        Self::new(vec![c])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
 
    #[test]
    fn append_arc() { 
        let mut t = Tng::empty();
        assert_eq!(t.ncomps(), 0);

        t.append_arc(Component::new(vec![0,1], false)); // [0-1]
        assert_eq!(t.ncomps(), 1);

        t.append_arc(Component::new(vec![2,3], false)); // [0-1] [2-3]
        assert_eq!(t.ncomps(), 2);

        t.append_arc(Component::new(vec![1,2], false)); // [0-1-2-3]
        assert_eq!(t.ncomps(), 1);

        t.append_arc(Component::new(vec![0,3], false)); // [-0-1-2-3-]
        assert_eq!(t.ncomps(), 1);
    }

    #[test]
    fn connect() { 
        let mut t0 = Tng::new(vec![
            Component::new(vec![0,1], false),
            Component::new(vec![2,3], false),
            Component::new(vec![10], true),
        ]);

        let t1 = Tng::new(vec![
            Component::new(vec![1,2], false),
            Component::new(vec![3,4], false),
            Component::new(vec![11], true),
        ]);

        t0.connect(t1);

        assert_eq!(t0, Tng::new(vec![
            Component::new(vec![0,1,2,3,4], false),
            Component::new(vec![10], true),
            Component::new(vec![11], true),
        ]));
    }

    #[test]
    fn deloop() { 
        let mut t = Tng::empty();
        assert_eq!(t.ncomps(), 0);
        assert_eq!(t.find_loop(), None);

        t.append_arc(Component::new(vec![0,1], false));
        assert_eq!(t.ncomps(), 1);
        assert_eq!(t.find_loop(), None);

        t.append_arc(Component::new(vec![2,3], false));
        assert_eq!(t.ncomps(), 2);
        assert_eq!(t.find_loop(), None);

        t.append_arc(Component::new(vec![2,3], false));
        assert_eq!(t.ncomps(), 2);
        assert_eq!(t.find_loop(), Some(1));

        t.remove_at(1);

        assert_eq!(t.ncomps(), 1);
        assert_eq!(t.find_loop(), None);
    }
}