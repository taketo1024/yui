use std::ops::Index;
use std::fmt;
use yui_utils::bitseq::{BitSeq, Bit};
use super::Resolution;

#[derive(Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct State(BitSeq);

impl State { 
    pub fn empty() -> Self { 
        State(BitSeq::empty())
    }

    pub fn zeros(l: usize) -> Self { 
        State(BitSeq::zeros(l))
    }

    pub fn ones(l: usize) -> Self { 
        State(BitSeq::ones(l))
    }

    pub fn is_empty(&self) -> bool { 
        self.0.is_empty()
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

    pub fn sub(&self, l: usize) -> State { 
        Self(self.0.sub(l))
    }

    pub fn is_sub(&self, other: &Self) -> bool { 
        self.0.is_sub(&other.0)
    }

    pub fn overwrite(&mut self, other: &Self) { 
        self.0.overwrite(&other.0)
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
    fn sub() { 
        let s = State::from_iter([1,0,0,1,1]);
        assert_eq!(s.sub(3), State::from_iter([1,0,0]));
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