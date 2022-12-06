use itertools::Itertools;
use num_traits::Pow;

use crate::links::links::{Link, State, Component};

use super::algebra::KhAlgGen;

#[derive(PartialEq, Eq, Hash, Debug)]
pub struct KhEnhState { 
    state: State,
    label: Vec<KhAlgGen>
}

impl KhEnhState {
    pub fn new(state: State, label: Vec<KhAlgGen>) -> KhEnhState { 
        KhEnhState { state, label }
    }
}

pub struct KhCubeVertex { 
    state: State,
    circles: Vec<Component>,
    generators: Vec<KhEnhState>
}

impl KhCubeVertex { 
    pub fn new(l: &Link, s: &State) -> Self {
        use super::algebra::KhAlgGen::{X, I};
        let circles = l.clone().resolve(&s).components();
        let r = circles.len();
        let generators = (0..2.pow(r)).map(|mut i| { 
            let labels = (0..r).map(|_j| { 
                let x = if i & 1 == 1 { I } else { X };
                i >>= 1;
                x
            }).collect();
            KhEnhState::new( s.clone(), labels )
        }).collect();

        KhCubeVertex { state: s.clone(), circles, generators }
    }
}

#[cfg(test)]
mod tests { 
    use crate::links::links::Resolution::*;
    use super::*;
    
    #[test]
    fn empty() { 
        let l = Link::empty();
        let s = State::empty();
        let v = KhCubeVertex::new(&l, &s);

        assert_eq!(v.state, s);
        assert_eq!(v.circles.len(), 0);
        assert_eq!(v.generators.len(), 1);
    }
    
    #[test]
    fn unknot() { 
        let l = Link::unknot();
        let s = State::empty();
        let v = KhCubeVertex::new(&l, &s);

        assert_eq!(v.state, s);
        assert_eq!(v.circles.len(), 1);
        assert_eq!(v.generators.len(), 2);
    }

    #[test]
    fn unlink_2() {
        let l = Link::from([[0, 1, 1, 0]]).resolve_at(0, Res1);
        let s = State::empty();
        let v = KhCubeVertex::new(&l, &s);

        assert_eq!(v.state, s);
        assert_eq!(v.circles.len(), 2);
        assert_eq!(v.generators.len(), 4);
    }
}