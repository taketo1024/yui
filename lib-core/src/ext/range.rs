use std::ops::{Add, Neg, Range, RangeInclusive, Sub};

pub trait RangeExt {
    type Idx; 
    fn mv(&self, l: Self::Idx, r: Self::Idx) -> Self;
    fn shift(&self, a: Self::Idx) -> Self;
    fn expand(&self, a: Self::Idx) -> Self;
}

impl<Idx> RangeExt for Range<Idx>
where Idx: Copy + Add<Output = Idx> + Sub<Output = Idx> + Neg<Output = Idx> {
    type Idx = Idx;

    fn mv(&self, l: Self::Idx, r: Self::Idx) -> Self {
        (self.start + l) .. (self.end + r)
    }

    fn shift(&self, a: Idx) -> Self {
        self.mv(a, a)
    }
    
    fn expand(&self, a: Self::Idx) -> Self {
        self.mv(-a, a)
    }    
}

impl<Idx> RangeExt for RangeInclusive<Idx>
where Idx: Copy + Add<Output = Idx> + Sub<Output = Idx> + Neg<Output = Idx> {
    type Idx = Idx;
    
    fn mv(&self, l: Self::Idx, r: Self::Idx) -> Self {
        (*self.start() + l) ..= (*self.end() + r)
    }

    fn shift(&self, a: Idx) -> Self {
        self.mv(a, a)
    }
    
    fn expand(&self, a: Self::Idx) -> Self {
        self.mv(-a, a)
    }    
}

#[cfg(test)]
mod tests {
    use crate::RangeExt;
 
    #[test]
    fn range() { 
        let r = -1 .. 3;
        assert_eq!(r.mv(2, 3), 1 .. 6);
        assert_eq!(r.shift(2), 1 .. 5);
        assert_eq!(r.expand(2), -3 .. 5);
    }

    #[test]
    fn range_incl() { 
        let r = -1 ..= 3;
        assert_eq!(r.mv(2, 3), 1 ..= 6);
        assert_eq!(r.shift(2), 1 ..= 5);
        assert_eq!(r.expand(2), -3 ..= 5);
    }
}