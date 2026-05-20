/// Division rounded to the nearest integer.
pub trait DivRound {
    fn div_round(&self, rhs: &Self) -> Self;
}