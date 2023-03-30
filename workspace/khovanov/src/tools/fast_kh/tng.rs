use std::fmt::Display;
use itertools::Itertools;
use yui_link::Component;

#[derive(Clone)]
pub struct Tng {
    comps: Vec<Component> // arc or circle
}

impl Tng { 
    pub fn new(comps: Vec<Component>) -> Self { 
        debug_assert!(comps.iter().all(|c| !c.is_empty()));
        Self { comps }
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

    pub fn comp(&self, i: usize) -> &Component { 
        &self.comps[i]
    }

    pub fn append_arc(&mut self, arc: Component) -> TngUpdate { 
        assert!(arc.is_arc());

        let n = self.ncomps();
        if let Some(i) = self.find_connectable(&arc, n) { 
            self.comps[i].connect(arc);

            // When both ends of `arc` are connectable:
            if let Some(j) = self.find_connectable(&self.comps[i], i) { 
                let arc_j = self.comps.remove(j);
                self.comps[i].connect(arc_j);
                TngUpdate(i, Some(j))
            } else { 
                TngUpdate(i, None)
            }
        } else { 
            self.comps.push(arc);
            TngUpdate(n, None)
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

    pub fn deloop(&mut self, i: usize) -> Component {
        assert!(self.comps[i].is_circle());
        self.comps.remove(i)
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
        let arcs = self.comps.iter().map(|a| 
            a.to_string()
        ).join(", ");
        write!(f, "T({})", arcs)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct TngUpdate(usize, Option<usize>);

impl TngUpdate { 
    pub(crate) fn new(index: usize, removed: Option<usize>) -> Self {
        Self(index, removed)
    }
    
    pub fn index(&self) -> usize { 
        self.0
    }

    pub fn removed(&self) -> Option<usize> { 
        self.1
    }

    pub fn reindex(&self, k: usize) -> usize { 
        let Some(j) = self.removed() else { return k };
        let i = self.index();

        if k < j { 
            k
        } else if k == j { 
            i
        } else {
            k - 1
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
 
    #[test]
    fn append() { 
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

        t.deloop(1);

        assert_eq!(t.ncomps(), 1);
        assert_eq!(t.find_loop(), None);
    }
}