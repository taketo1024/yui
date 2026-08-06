//! [`IntoDigits`]: decompose an integer into its base-10 digits.

/// Decompose a non-negative integer into its base-10 digits.
///
/// See: <https://en.wikipedia.org/wiki/Positional_notation>
pub trait IntoDigits: Sized {
    type Digit;

    /// Iterator over digits, most-significant first.
    fn into_digits(self) -> impl Iterator<Item = Self::Digit> {
        let mut v: Vec<_> = self.into_rev_digits().collect();
        v.reverse();
        v.into_iter()
    }

    /// Iterator over digits, least-significant first.
    fn into_rev_digits(self) -> impl Iterator<Item = Self::Digit> {
        let mut v: Vec<_> = self.into_digits().collect();
        v.reverse();
        v.into_iter()
    }
}

macro_rules! impl_into_digits {
    ($t: ty, $d: ty) => {
        impl IntoDigits for $t {
            type Digit = $d;
            fn into_rev_digits(self) -> impl Iterator<Item = $d> {
                std::iter::successors(
                    Some(self),
                    |&n| (n >= 10).then(|| n / 10),
                ).map(|n| (n % 10) as $d)
            }
        }
    };
}

impl_into_digits!(usize, u8);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn usize() {
        let a = 123456789;
        assert_eq!(a.into_digits().collect::<Vec<_>>(), vec![1,2,3,4,5,6,7,8,9]);
        assert_eq!(a.into_rev_digits().collect::<Vec<_>>(), vec![9,8,7,6,5,4,3,2,1]);
    }

    #[test]
    fn zero() {
        assert_eq!(0usize.into_digits().collect::<Vec<_>>(), vec![0]);
        assert_eq!(0usize.into_rev_digits().collect::<Vec<_>>(), vec![0]);
    }
}
