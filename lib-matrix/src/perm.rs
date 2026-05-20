//! Permutation of `0..n`, used to track row/column reorderings during
//! matrix decompositions (PLUQ, pivot search, Schur complement).
//!
//! Stored as a `Vec<usize>` where `perm[i]` is the image of `i`.

use std::ops::{Mul, MulAssign};
use auto_impl_ops::auto_ops;

/// A permutation of `0..n`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Perm {
    data: Vec<usize>,
}

impl Perm {
    /// Create from the image vector. `data[i]` is the image of `i`.
    /// Panics (debug) if `data` is not a valid permutation of `0..data.len()`.
    pub fn new(data: Vec<usize>) -> Self {
        debug_assert!(is_valid_perm(&data), "not a valid permutation: {:?}", data);
        Self { data }
    }

    /// The identity permutation on `0..n`.
    pub fn id(n: usize) -> Self {
        Self { data: (0..n).collect() }
    }

    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Alias for [`len`](Self::len), matching `sprs`'s naming.
    pub fn dim(&self) -> usize {
        self.data.len()
    }

    /// The image of `i`.
    pub fn at(&self, i: usize) -> usize {
        self.data[i]
    }

    /// `true` iff this is the identity permutation.
    pub fn is_id(&self) -> bool {
        self.data.iter().enumerate().all(|(i, &x)| i == x)
    }

    /// The inverse permutation.
    pub fn inv(&self) -> Self {
        let n = self.data.len();
        let mut inv = vec![0; n];
        for (i, &j) in self.data.iter().enumerate() {
            inv[j] = i;
        }
        Self { data: inv }
    }

    /// Underlying image slice (`raw()[i]` is the image of `i`).
    pub fn raw(&self) -> &[usize] {
        &self.data
    }

    /// Permutation `p` of `0..n` that sends each index in `prefix` to a
    /// position `0, 1, 2, ...` (in the order given), with the remaining
    /// indices filling positions in sorted order.
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
        Perm::new(rhs.data.iter().map(|&i| self.at(i)).collect())
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
        assert!(!Perm::new(vec![1, 0]).is_id());
    }

    #[test]
    fn accessors() {
        let v = vec![2, 0, 1, 3];
        let p = Perm::new(v.clone());
        assert_eq!(p.len(), 4);
        assert_eq!(p.dim(), p.len());
        assert_eq!(p.raw(), v.as_slice());
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

    // --- auto_ops-derived variants ---

    // --- pull_and_fill ---

    #[test]
    fn pull_and_fill_basic() {
        // n=5, prefix=[3,1] → sends 3→0, 1→1, others fill sorted: 0→2, 2→3, 4→4.
        let p = Perm::pull_and_fill(5, [3, 1]);
        assert_eq!(p.raw(), &[2, 1, 3, 0, 4]);
        assert_eq!(p.at(3), 0);
        assert_eq!(p.at(1), 1);
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
