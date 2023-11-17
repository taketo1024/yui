// Makes a k-repetitive combination of indices in 0..n 
// from an n-combination of indices in 0..(n + k). 
//
// e.g. (n, k) = (4, 7)
//
//  [1,5,7,10] = *|***|*|**|
//   -> [(0,1), (1,3), (2,1), (3,2)]
//
// The last element of the input represents the right-end wall,
// which gives (n + k - 1).

pub fn rep_comb(list: &[usize]) -> Vec<(usize, usize)> { 
    let mut i = 0;
    let mut s = 0; // starting index for i.

    list.iter().filter_map(|&p| { 
        let e = if p > s {
            Some((i, p - s))
        } else { 
            None
        };
        i += 1;
        s = p + 1;
        e
    }).collect()
}

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn sample() { 
        let a = rep_comb(&[1,5,7,10]);
        assert_eq!(a, vec![(0,1), (1,3), (2,1), (3,2)]);
    }

    #[test]
    fn empty() { 
        let a = rep_comb(&[]);
        assert!(a.is_empty());
    }

    #[test]
    fn only_walls() { 
        let a = rep_comb(&[0,1,2,3]);
        assert!(a.is_empty());
    }

    #[test]
    fn large_value() { 
        let a = rep_comb(&[usize::MAX-1]);
        assert_eq!(a, vec![(0, usize::MAX-1)]);
    }
}