use std::collections::HashSet;
use std::iter::zip;
use itertools::Itertools;
use yui_link::{Edge, Link, Path};

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
    fn colored_seifert_circles(&self, base: Edge) -> Vec<(Path, Color)>;
}

impl LinkExt for Link { 
    fn colored_seifert_circles(&self, base: Edge) -> Vec<(Path, Color)> {
        assert_eq!(self.components().len(), 1, "Only knots are supported.");

        let circles = self.seifert_circles();
        let n = circles.len();
    
        let mut colors = vec![Color::A; n];
        let mut queue = vec![];
        let mut remain: HashSet<_> = (0..n).collect();
    
        let i = circles.iter().find_position(|c| 
            c.edges().contains(&base)
        ).unwrap().0;
    
        queue.push(i);
        colors[i] = Color::A;
    
        while !queue.is_empty() { 
            let i1 = queue.remove(0);
            let c1 = &circles[i1];
    
            let adjs = remain.iter().filter_map(|&i2| {
                let c2 = &circles[i2];
                if c1.is_adj(c2, self) { Some(i2) } else { None }
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