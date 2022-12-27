use std::vec::IntoIter;
use itertools::Itertools;
use derive_more::{Display, Add, Sub};

#[derive(Clone, Copy, PartialEq, Eq, Hash, Display, Add, Sub)]
#[display(fmt = "({}, {})", _0, _1)]
pub struct Idx2(pub isize, pub isize);

impl Idx2 { 
    pub fn iterate(from: Idx2, to: Idx2, step:(usize, usize)) -> Idx2Range {
        (from.1 ..= to.1).step_by(step.1).flat_map(|j| { 
            (from.0 ..= to.0).step_by(step.0).map(move |i| Idx2(i, j))
        }).collect_vec().into_iter()
    }

    pub fn as_tuple(&self) -> (isize, isize) {
        (self.0, self.1)
    }
}

impl From<[isize; 2]> for Idx2 { 
    fn from(idx: [isize; 2]) -> Self {
        Self(idx[0], idx[1])
    }
}

pub type Idx2Range = IntoIter<Idx2>;