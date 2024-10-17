use std::ops::{Add, Range, RangeInclusive, Sub};

pub trait RangeExt {
    type Idx; 
    fn shift(&self, a: Self::Idx) -> Self;
    fn extend(&self, a: Self::Idx) -> Self;
}

impl<Idx> RangeExt for Range<Idx>
where Idx: Copy + Add<Output = Idx> + Sub<Output = Idx> {
    type Idx = Idx;
    fn shift(&self, a: Idx) -> Self {
        (self.start + a) .. (self.end + a)
    }
    
    fn extend(&self, a: Self::Idx) -> Self {
        (self.start - a) .. (self.end + a)
    }
}

impl<Idx> RangeExt for RangeInclusive<Idx>
where Idx: Copy + Add<Output = Idx> + Sub<Output = Idx> {
    type Idx = Idx;
    fn shift(&self, a: Idx) -> Self {
        (*self.start() + a) ..= (*self.end() + a)
    }
    
    fn extend(&self, a: Self::Idx) -> Self {
        (*self.start() - a) ..= (*self.end() + a)
    }
}

#[cfg(test)]
mod tests {
    use crate::RangeExt;
 
    #[test]
    fn range() { 
        let r = -1 .. 3;
        assert_eq!(r.shift(2), 1 .. 5);
        assert_eq!(r.extend(2), -3 .. 5);
    }

    #[test]
    fn range_incl() { 
        let r = -1 ..= 3;
        assert_eq!(r.shift(2), 1 ..= 5);
        assert_eq!(r.extend(2), -3 ..= 5);
    }
}