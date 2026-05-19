use is_even::IsEven;
use num_traits::One;

pub trait PowMod2<Rhs> 
where Self: One, Rhs: IsEven {
    fn pow_mod2(self, rhs: Rhs) -> Self { 
        if rhs.is_even() { Self::one() } else { self }
    }
}

impl<T, Rhs> PowMod2<Rhs> for T
where T: One, Rhs: IsEven {}

