//! Permutation of `0..n`, used to track row/column reorderings during
//! matrix decompositions (PLUQ, pivot search, Schur complement).
//!
//! Internally, the identity permutation is represented by just its size,
//! making `Perm::id(n)` zero-cost.

use std::ops::{Mul, MulAssign};
use auto_impl_ops::auto_ops;
use either::Either;

/// A permutation of `0..n`.
///
/// Stored as `Either<usize, Vec<usize>>`: `Left(n)` is the identity
/// permutation on `0..n`, and `Right(v)` is an explicit image vector
/// where `v[i]` is the image of `i`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Perm {
    data: Either<usize, Vec<usize>>,
}

impl Perm {
    /// Create from the image vector. `data[i]` is the image of `i`.
    /// Panics (debug) if `data` is not a valid permutation of `0..data.len()`.
    pub(crate) fn new(data: Vec<usize>) -> Self {
        debug_assert!(is_valid_perm(&data), "not a valid permutation: {:?}", data);
        Self { data: Either::Right(data) }
    }

    /// Create from an iterator of images. The `i`-th item is the image of `i`.
    /// Panics (debug) if the resulting sequence is not a valid permutation.
    pub fn from_indices<I>(images: I) -> Self
    where I: IntoIterator<Item = usize> {
        Self::new(images.into_iter().collect())
    }

    /// The identity permutation on `0..n`. Zero-cost.
    pub fn id(n: usize) -> Self {
        Self { data: Either::Left(n) }
    }

    pub fn len(&self) -> usize {
        match &self.data {
            Either::Left(n) => *n,
            Either::Right(v) => v.len(),
        }
    }

    /// Alias for [`len`](Self::len), matching `sprs`'s naming.
    pub fn dim(&self) -> usize {
        self.len()
    }

    /// The image of `i`.
    pub fn at(&self, i: usize) -> usize {
        match &self.data {
            Either::Left(_) => i,
            Either::Right(v) => v[i],
        }
    }

    /// `true` iff this is the identity permutation.
    pub fn is_id(&self) -> bool {
        match &self.data {
            Either::Left(_) => true,
            Either::Right(v) => v.iter().enumerate().all(|(i, &x)| i == x),
        }
    }

    /// The inverse permutation.
    pub fn inv(&self) -> Self {
        match &self.data {
            Either::Left(n) => Self::id(*n),
            Either::Right(v) => {
                let mut inv = vec![0; v.len()];
                for (i, &j) in v.iter().enumerate() {
                    inv[j] = i;
                }
                Self { data: Either::Right(inv) }
            }
        }
    }

    /// Permutation `p` of `0..n` that sends each index in `prefix` to a
    /// position `0, 1, 2, ...` (in the order given), with the remaining
    /// indices filling positions in sorted order.
    /// If `prefix` is empty, returns the identity (zero-cost).
    pub fn pull_and_fill<I>(n: usize, prefix: I) -> Self
    where I: IntoIterator<Item = usize> {
        let mut data = vec![0usize; n];
        let mut taken = vec![false; n];
        let mut k = 0;
        for i in prefix {
            debug_assert!(i < n, "index {i} out of range 0..{n}");
            debug_assert!(!taken[i], "duplicate index {i} in prefix");
            data[i] = k;
            taken[i] = true;
            k += 1;
        }
        if k == 0 {
            return Self::id(n);
        }
        let mut pos = k;
        for j in 0..n {
            if !taken[j] {
                data[j] = pos;
                pos += 1;
            }
        }
        Self::new(data)
    }
}

/// Composition `(p * q)(i) = p(q(i))` (right-to-left, math convention).
#[auto_ops]
impl<'a, 'b> Mul<&'b Perm> for &'a Perm {
    type Output = Perm;
    fn mul(self, rhs: &'b Perm) -> Perm {
        assert_eq!(self.len(), rhs.len(), "permutations must have the same length");
        match (&self.data, &rhs.data) {
            (Either::Left(_), _) => rhs.clone(),
            (_, Either::Left(_)) => self.clone(),
            (Either::Right(_), Either::Right(qv)) => {
                Perm::new(qv.iter().map(|&i| self.at(i)).collect())
            }
        }
    }
}

fn is_valid_perm(v: &[usize]) -> bool {
    let n = v.len();
    let mut seen = vec![false; n];
    for &x in v {
        if x >= n || seen[x] {
            return false;
        }
        seen[x] = true;
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    // --- constructors ---

    #[test]
    fn id() {
        let p = Perm::id(4);
        assert!(p.is_id());
        assert_eq!(p.len(), 4);
        for i in 0..4 {
            assert_eq!(p.at(i), i);
        }
    }

    #[test]
    fn new_accepts_valid() {
        let p = Perm::new(vec![2, 0, 1]);
        assert_eq!(p.at(0), 2);
        assert_eq!(p.at(1), 0);
        assert_eq!(p.at(2), 1);
    }

    #[test]
    #[should_panic]
    fn new_rejects_duplicate() {
        let _ = Perm::new(vec![0, 0, 1]);
    }

    #[test]
    #[should_panic]
    fn new_rejects_out_of_range() {
        let _ = Perm::new(vec![0, 1, 5]);
    }

    // --- accessors ---

    #[test]
    fn is_id_true_false() {
        assert!(Perm::id(3).is_id());
        // identity image vector via `new` is still detected as identity
        assert!(Perm::new(vec![0, 1, 2]).is_id());
        assert!(!Perm::new(vec![1, 0]).is_id());
    }

    #[test]
    fn accessors() {
        let v = vec![2, 0, 1, 3];
        let p = Perm::new(v.clone());
        assert_eq!(p.len(), 4);
        assert_eq!(p.dim(), p.len());
        for (i, &expected) in v.iter().enumerate() {
            assert_eq!(p.at(i), expected);
        }
    }

    // --- inverse ---

    #[test]
    fn inv_of_id() {
        let id = Perm::id(5);
        assert_eq!(id.inv(), id);
    }

    #[test]
    fn inv_of_inv() {
        let p = Perm::new(vec![2, 0, 3, 1]);
        assert_eq!(p.inv().inv(), p);
    }

    #[test]
    fn inv_roundtrip() {
        let p = Perm::new(vec![2, 0, 3, 1]);
        let pi = p.inv();
        for i in 0..p.len() {
            assert_eq!(pi.at(p.at(i)), i);
            assert_eq!(p.at(pi.at(i)), i);
        }
    }

    // --- composition ---

    #[test]
    fn mul_formula() {
        // (p * q)(i) = p(q(i))
        let p = Perm::new(vec![2, 0, 3, 1]);
        let q = Perm::new(vec![1, 2, 0, 3]);
        let pq = &p * &q;
        for i in 0..4 {
            assert_eq!(pq.at(i), p.at(q.at(i)));
        }
    }

    #[test]
    fn mul_with_id() {
        let p = Perm::new(vec![2, 0, 3, 1]);
        let id = Perm::id(4);
        assert_eq!(&p * &id, p);
        assert_eq!(&id * &p, p);
    }

    #[test]
    fn mul_id_with_id() {
        let id = Perm::id(4);
        assert_eq!(&id * &id, id);
    }

    #[test]
    fn mul_by_inverse_is_id() {
        let p = Perm::new(vec![2, 0, 3, 1]);
        assert!((&p * &p.inv()).is_id());
        assert!((&p.inv() * &p).is_id());
    }

    #[test]
    fn mul_associative() {
        let p = Perm::new(vec![2, 0, 3, 1]);
        let q = Perm::new(vec![1, 3, 0, 2]);
        let r = Perm::new(vec![3, 1, 2, 0]);
        assert_eq!(&(&p * &q) * &r, &p * &(&q * &r));
    }

    #[test]
    #[should_panic]
    fn mul_panics_on_dim_mismatch() {
        let p = Perm::id(3);
        let q = Perm::id(4);
        let _ = &p * &q;
    }

    // --- pull_and_fill ---

    #[test]
    fn pull_and_fill_basic() {
        // n=5, prefix=[3,1] → sends 3→0, 1→1, others fill sorted: 0→2, 2→3, 4→4.
        let p = Perm::pull_and_fill(5, [3, 1]);
        let expected = [2, 1, 3, 0, 4];
        for (i, &x) in expected.iter().enumerate() {
            assert_eq!(p.at(i), x);
        }
    }

    #[test]
    fn pull_and_fill_empty_prefix() {
        let p = Perm::pull_and_fill(4, std::iter::empty());
        assert!(p.is_id());
    }

    #[test]
    fn pull_and_fill_full_prefix() {
        // Specifying every index reduces to: p(prefix[k]) = k.
        let p = Perm::pull_and_fill(4, [2, 0, 3, 1]);
        assert_eq!(p.at(2), 0);
        assert_eq!(p.at(0), 1);
        assert_eq!(p.at(3), 2);
        assert_eq!(p.at(1), 3);
    }

    // --- auto_ops-derived variants ---

    #[test]
    fn mul_all_ref_variants() {
        let p = Perm::new(vec![2, 0, 3, 1]);
        let q = Perm::new(vec![1, 2, 0, 3]);
        let expected = &p * &q;

        assert_eq!(p.clone() * q.clone(), expected);
        assert_eq!(p.clone() * &q, expected);
        assert_eq!(&p * q.clone(), expected);
        assert_eq!(&p * &q, expected);
    }
}
