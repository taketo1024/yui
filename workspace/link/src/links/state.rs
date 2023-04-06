use std::ops::Index;
use std::fmt;
use yui_utils::bitseq::{BitSeq, Bit};
use super::Resolution;

#[derive(Debug, Clone, Default, PartialEq, Eq, Hash, Ord)]
pub struct State(BitSeq);

impl State { 
    pub fn empty() -> Self { 
        State(BitSeq::empty())
    }

    pub fn zeros(l: usize) -> Self { 
        assert!(l <= BitSeq::MAX_LEN);
        State(BitSeq::new(0, l))
    }

    pub fn ones(l: usize) -> Self { 
        assert!(l <= BitSeq::MAX_LEN);
        let b = (0..l).fold(0, |b, _| b << 1 | 1);
        State(BitSeq::new(b, l))
    }

    pub fn len(&self) -> usize { 
        self.0.len()
    }

    pub fn weight(&self) -> usize { 
        self.0.weight()
    }

    pub fn iter(&self) -> impl Iterator<Item = Resolution> { 
        self.0.iter().map(|b| 
            if b.is_zero() { 
                Resolution::Res0
            } else { 
                Resolution::Res1
            }
        )
    }
    
    pub fn push(&mut self, r: Resolution) { 
        if r.is_zero() { 
            self.0.push(Bit::Bit0)
        } else { 
            self.0.push(Bit::Bit1)
        }
    }

    pub fn append(&mut self, other: Self) {
        self.0.append(other.0);
    }

    pub fn targets(&self) -> Vec<State> { 
        let n = self.len();
        (0..n).filter(|&i| self[i].is_zero() ).map(|i| { 
            let mut t = self.clone();
            t.0.set_1(i);
            t
        }).collect()
    }

    pub fn generate(len: usize) -> Vec<State> { 
        BitSeq::generate(len).into_iter().map(|b| 
            Self(b)
        ).collect()
    }
}

impl From<BitSeq> for State {
    fn from(b: BitSeq) -> Self {
        Self(b)
    }
}

impl<T> FromIterator<T> for State
where Bit: From<T> {
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        Self(BitSeq::from_iter(iter))
    }
}

impl Index<usize> for State { 
    type Output = Resolution;
    fn index(&self, index: usize) -> &Resolution {
        if self.0[index].is_zero() { 
            &Resolution::Res0
        } else { 
            &Resolution::Res1
        }
    }
}

impl fmt::Display for State {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.0.fmt(f)
    }
}

impl PartialOrd for State {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        let cmp = self.len().cmp(&other.len()).then_with( ||
            self.weight().cmp(&other.weight())
        ).then_with( ||
            self.0.cmp(&other.0).reverse()
        );
        Some(cmp)
    }
}

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn state_targets() {
        let s = State::from_iter([0, 0, 0]);
        assert_eq!(s.targets(), vec![
            State::from_iter([1, 0, 0]),
            State::from_iter([0, 1, 0]),
            State::from_iter([0, 0, 1]),
        ]);

        let s = State::from_iter([0, 1, 0]);
        assert_eq!(s.targets(), vec![
            State::from_iter([1, 1, 0]),
            State::from_iter([0, 1, 1]),
        ]);

        let s = State::from_iter([1, 1, 0]);
        assert_eq!(s.targets(), vec![
            State::from_iter([1, 1, 1]),
        ]);

        let s = State::from_iter([1, 1, 1]);
        assert_eq!(s.targets(), vec![]);
    }

    #[test]
    fn order() { 
        let mut ss = State::generate(3);
        ss.sort();
        
        assert_eq!(ss, 
            [[0, 0, 0],[1, 0, 0],[0, 1, 0],[0, 0, 1],[1, 1, 0],[1, 0, 1],[0, 1, 1],[1, 1, 1]]
            .map(|v| State::from_iter(v))
        );
    }
}