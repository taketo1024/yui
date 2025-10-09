use std::ops::Index;
use std::fmt;

#[derive(Clone, PartialEq, Eq)]
pub struct Partition {
    parts: Vec<usize>,
}

impl Partition {
    /// Creates a new partition from a vector of non-negative integers.
    /// The parts are expected to be in non-increasing order.
    pub fn new(parts: Vec<usize>) -> Self {
        if !parts.windows(2).all(|w| w[0] >= w[1]) {
            panic!("Partition parts must be in non-increasing order");
        }
        Partition { parts }
    }

    /// Returns a reference to the parts of the partition.
    pub fn parts(&self) -> &[usize] {
        &self.parts
    }

    /// Returns the sum of the parts (the integer being partitioned).
    pub fn sum(&self) -> usize {
        self.parts.iter().sum()
    }

    /// Returns the number of parts in the partition.
    pub fn len(&self) -> usize {
        self.parts.len()
    }

    /// Returns true if the partition has no parts.
    pub fn is_empty(&self) -> bool {
        self.parts.is_empty()
    }

    /// Returns true if the Young diagram of self contains that of partition `p`.
    pub fn contains(&self, p: &Partition) -> bool {
        (0..p.len()).all(|i| self[i] >= p[i])
    }

    /// Generates all partitions of a non-negative integer `n`.
    /// Returns an iterator over all partitions of a non-negative integer `n`.
    pub fn all_partitions(n: usize) -> PartitionIter {
        PartitionIter::new(n)
    }
}

impl<I> From<I> for Partition
where
    I: IntoIterator<Item = usize>,
{
    fn from(iter: I) -> Self {
        let parts = iter.into_iter().collect();
        Partition::new(parts)
    }
}

impl FromIterator<usize> for Partition {
    fn from_iter<T: IntoIterator<Item = usize>>(iter: T) -> Self {
        Partition::from(iter)
    }
}

impl Index<usize> for Partition {
    type Output = usize;

    fn index(&self, index: usize) -> &Self::Output {
        self.parts.get(index).unwrap_or(&0)
    }
}

impl fmt::Display for Partition {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "(")?;
        for (i, part) in self.parts.iter().enumerate() {
            if i > 0 {
                write!(f, ", ")?;
            }
            write!(f, "{}", part)?;
        }
        write!(f, ")")
    }
}

impl fmt::Debug for Partition {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_tuple("Partition")
            .field(&self.parts)
            .finish()
    }
}

impl Partition {
    /// Returns the Young diagram of the partition as a String.
    /// Each row is represented by a line of '■' characters.
    pub fn young_diagram(&self) -> String {
        self.parts
            .iter()
            .map(|&n| "⬜︎".repeat(n))
            .collect::<Vec<_>>()
            .join("\n")
    }
}

/// Iterator over all partitions of a non-negative integer n.
pub struct PartitionIter {
    a: Vec<usize>,
    k: usize,
    done: bool,
}

impl PartitionIter {
    pub fn new(n: usize) -> Self {
        let mut a = vec![0; n + 1];
        a[0] = n;
        let k = if n > 0 { 1 } else { 0 };
        
        PartitionIter {
            a,
            k,
            done: false,
        }
    }
}

impl Iterator for PartitionIter {
    type Item = Partition;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        // First call: return the initial partition
        let k = self.k;
        let result = Partition::new(self.a[..k].to_vec());

        // Generate the next partition using the standard algorithm
        // (see e.g. https://jeromekelleher.net/generating-integer-partitions.html)

        // Find the lowest index whose part is greater than 1
        if let Some(i) = (0..k).rev().find(|&i| self.a[i] > 1) { 
            self.a[i] -= 1;

            let a = self.a[i];
            let mut rem_val = k - i; // (1's between (i + 1)..k) + 1

            for j in i + 1 .. { 
                if rem_val > a {
                    self.a[j] = a;
                    rem_val -= a;
                } else { 
                    self.a[j] = rem_val;
                    self.k = j + 1;
                    break;
                }
            }
        } else { 
            self.done = true
        }

        Some(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_partition_valid() {
        let p = Partition::new(vec![4, 2, 1]);
        assert_eq!(p.parts(), &[4, 2, 1]);
    }

    #[test]
    #[should_panic(expected = "Partition parts must be in non-increasing order")]
    fn test_new_partition_invalid_order() {
        Partition::new(vec![2, 4, 1]);
    }

    #[test]
    fn test_from() {
        let v = vec![3, 2, 1];
        let p = Partition::from(v);
        assert_eq!(p.parts(), &[3, 2, 1]);
    }

    #[test]
    fn test_from_iterator() {
        let v = vec![3, 2, 1];
        let p: Partition = v.iter().cloned().collect::<Vec<_>>().into();
        assert_eq!(p.parts(), &[3, 2, 1]);
    }

    #[test]
    fn test_sum() {
        let p = Partition::new(vec![3, 2, 1]);
        assert_eq!(p.sum(), 6);
    }

    #[test]
    fn test_len_and_is_empty() {
        let p = Partition::new(vec![2, 2]);
        assert_eq!(p.len(), 2);
        assert!(!p.is_empty());

        let empty = Partition::new(vec![]);
        assert_eq!(empty.len(), 0);
        assert!(empty.is_empty());
    }

    #[test]
    fn test_index() {
        let p = Partition::new(vec![5, 3, 1]);
        assert_eq!(p[0], 5);
        assert_eq!(p[1], 3);
        assert_eq!(p[2], 1);
        assert_eq!(p[3], 0); // out of bounds returns 0
    }

    #[test]
    fn test_display() {
        let p = Partition::new(vec![3, 2, 1]);
        assert_eq!(format!("{}", p), "(3, 2, 1)");
    }

    #[test]
    fn test_debug() {
        let p = Partition::new(vec![2, 1]);
        assert_eq!(format!("{:?}", p), "Partition([2, 1])");
    }

    #[test]
    fn test_all_partitions() {
        let parts = Partition::all_partitions(4);
        let expected: Vec<Vec<usize>> = vec![
            vec![4],
            vec![3, 1],
            vec![2, 2],
            vec![2, 1, 1],
            vec![1, 1, 1, 1],
        ];
        let actual: Vec<Vec<usize>> = parts.into_iter().map(|p| p.parts().to_vec()).collect();
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_all_partitions_zero() {
        let mut iter = Partition::all_partitions(0);
        let p = iter.next();
        assert_eq!(p, Some(Partition::from([])));
        assert!(iter.next().is_none());
    }

    #[test]
    fn test_all_partitions_one() {
        let mut iter = Partition::all_partitions(1);
        let p = iter.next();
        assert_eq!(p.unwrap().parts(), &[1]);
        assert!(iter.next().is_none());
    }

    #[test]
    fn test_contains() {
        let p = Partition::new(vec![5, 3, 2]);
        let q = Partition::new(vec![4, 2, 1]);
        let r = Partition::new(vec![5, 3, 2]);
        let s = Partition::new(vec![6, 3, 2]);
        let t = Partition::new(vec![5, 4, 2]);
        let u = Partition::new(vec![5, 3, 2, 1]);
        let v = Partition::new(vec![]);

        // p contains q (all parts of p >= q)
        assert!(p.contains(&q));
        // p contains itself
        assert!(p.contains(&r));
        // p does not contain s (5 < 6)
        assert!(!p.contains(&s));
        // p does not contain t (3 < 4)
        assert!(!p.contains(&t));
        // p does not contain u (u has more parts)
        assert!(!p.contains(&u));
        // Any partition contains the empty partition
        assert!(p.contains(&v));
        // The empty partition contains only itself
        assert!(v.contains(&v));
        assert!(!v.contains(&p));
    }

    #[test]
    fn test_young_diagram() {
        let p = Partition::new(vec![3, 2, 1]);
        let expected = "⬜︎⬜︎⬜︎\n⬜︎⬜︎\n⬜︎";
        assert_eq!(p.young_diagram(), expected);
    }
}
