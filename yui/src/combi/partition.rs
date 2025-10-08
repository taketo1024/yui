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

    /// Generates all partitions of a non-negative integer `n`.
    pub fn all_partitions(n: usize) -> Vec<Partition> {
        fn helper(n: usize, max: usize, prefix: &mut Vec<usize>, result: &mut Vec<Partition>) {
            if n == 0 {
                result.push(Partition::new(prefix.clone()));
                return;
            }
            for i in (1..=std::cmp::min(n, max)).rev() {
                prefix.push(i);
                helper(n - i, i, prefix, result);
                prefix.pop();
            }
        }

        let mut result = Vec::new();
        let mut prefix = Vec::new();
        helper(n, n, &mut prefix, &mut result);
        result
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
            .map(|&n| "◻︎".repeat(n))
            .collect::<Vec<_>>()
            .join("\n")
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
    fn test_young_diagram() {
        let p = Partition::new(vec![3, 2, 1]);
        let expected = "◻︎◻︎◻︎\n◻︎◻︎\n◻︎";
        assert_eq!(p.young_diagram(), expected);

        for p in Partition::all_partitions(5) { 
            println!("{p}\n{}\n", p.young_diagram())
        }
    }
}
