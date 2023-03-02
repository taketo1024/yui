use std::ops::Index;
use std::fmt;
use itertools::{join, Itertools};
use super::Resolution;
use Resolution::{Res0, Res1};

#[derive(Debug, Clone, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct State { 
    values: Vec<Resolution>
}

impl State { 
    pub fn empty() -> Self { 
        State { values: vec![] }
    }

    pub fn from_iter<I>(iter: I) -> Self
    where I: Iterator<Item = u8> {
        let values = iter.map(|a| Resolution::from(a)).collect();
        State{ values }
    }

    pub fn from_bseq(mut bseq: usize, length: usize) -> Self {
        let seq = (0..length)
            .map(|_| {
                let a = (bseq & 1) as u8;
                bseq >>= 1;
                Resolution::from(a)
            })
            .rev()
            .collect();
        State{ values: seq }
    }

    pub fn values(&self) -> &Vec<Resolution> {
        &self.values
    }

    pub fn weight(&self) -> usize { 
        self.values.iter().filter(|r| !r.is_zero()).count()
    }
    
    pub fn len(&self) -> usize { 
        self.values.len()
    }

    pub fn targets(&self) -> Vec<State> { 
        let n = self.len();
        (0..n).filter(|&i| self[i] == Res0 ).map(|i| { 
            let mut t = self.clone();
            t.values[i] = Res1;
            t
        }).collect_vec()
    }

    pub fn append(&mut self, mut other: Self) {
        self.values.append(&mut other.values);
    }
}

impl<const N: usize> From<[u8; N]> for State {
    fn from(values: [u8; N]) -> Self {
        Self::from_iter(values.into_iter())
    }
}

impl From<Vec<u8>> for State { 
    fn from(v: Vec<u8>) -> Self {
        Self::from_iter(v.into_iter())
    }
}

impl Index<usize> for State { 
    type Output = Resolution;
    fn index(&self, index: usize) -> &Resolution {
        &self.values[index]
    }
}

impl fmt::Display for State {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "[{}]", join(self.values.iter(), ", "))
    }
}

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn state_targets() {
        let s = State::from(vec![0, 0, 0]);
        assert_eq!(s.targets(), vec![
            State::from(vec![1, 0, 0]),
            State::from(vec![0, 1, 0]),
            State::from(vec![0, 0, 1]),
        ]);

        let s = State::from(vec![0, 1, 0]);
        assert_eq!(s.targets(), vec![
            State::from(vec![1, 1, 0]),
            State::from(vec![0, 1, 1]),
        ]);

        let s = State::from(vec![1, 1, 0]);
        assert_eq!(s.targets(), vec![
            State::from(vec![1, 1, 1]),
        ]);

        let s = State::from(vec![1, 1, 1]);
        assert_eq!(s.targets(), vec![]);
    }
}