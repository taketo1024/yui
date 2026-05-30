use std::collections::HashSet;
use std::iter::zip;
use itertools::Itertools;
use yui_link::{Link, Path};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Color { A, B }

impl Color { 
    pub fn is_a(&self) -> bool { 
        self == &Color::A
    }

    pub fn other(&self) -> Self { 
        match self { 
            Color::A => Color::B,
            Color::B => Color::A
        }
    }
}

pub trait LinkExt { 
    fn colored_seifert_circles(&self) -> Vec<(Path, Color)>;
}

impl LinkExt for Link { 
    fn colored_seifert_circles(&self) -> Vec<(Path, Color)> {
        assert!(self.is_knot(), "Only knots are supported.");
        assert!(self.base_pt().is_some());

        let circles = self.seifert_circles();
        let base_pt = self.base_pt().unwrap();
        let n = circles.len();
    
        let mut colors = vec![Color::A; n];
        let mut queue = vec![];
        let mut remain: HashSet<_> = (0..n).collect();
    
        let i = circles.iter().find_position(|c| 
            c.edges().contains(&base_pt)
        ).unwrap().0;
    
        queue.push(i);
        colors[i] = Color::A;
    
        while !queue.is_empty() { 
            let i1 = queue.remove(0);
            let c1 = &circles[i1];
    
            let adjs = remain.iter().filter_map(|&i2| {
                let c2 = &circles[i2];
                if is_adj(c1, c2, self) { Some(i2) } else { None }
            }).collect_vec();
            
            for i2 in adjs {
                remain.remove(&i2);
                queue.push(i2);
                colors[i2] = colors[i1].other();
            };
        }
    
        assert!(queue.is_empty());
    
        zip(circles, colors).collect()
    }
}

/// Two paths are adjacent if some crossing of `link` touches an edge on
/// `self` and a *different* edge on `other`.
fn is_adj(p1: &Path, p2: &Path, link: &Link) -> bool {
    for x in link.nodes() {
        if !x.edges().iter().any(|e| p1.contains(*e)) {
            continue
        }

        let Some(&e) = x.edges().iter().find(|e| !p1.contains(**e)) else {
            continue
        };

        if p2.contains(e) {
            return true
        }
    }

    false
}