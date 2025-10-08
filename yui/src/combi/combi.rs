use dinglebit_combinatorics::Combination;

pub fn combi(n: usize, r: usize) -> impl Iterator<Item = Vec<usize>> {
    Combination::new(n, r).into_iter()
}

pub fn rep_combi_count(n: usize, r: usize) -> impl Iterator<Item = Vec<(usize, usize)>> {
    // e.g. when (n, k) = (4, 7),
    // part: [1,5,7] = *|***|*|**
    //                  1   5 7
    //    => [(0,1), (1,3), (2,1), (3,2)] 

    #[allow(non_snake_case)]
    let (N, k) = if n > 0 { 
        (n + r - 1, n - 1) 
    } else if r > 0 { 
        (0, 1)
    } else { 
        (0, 0)
    };

    combi(N, k).map(move |mut part| { 
        part.push(N); // right-end wall
        part.into_iter().fold((0, 0, vec![]), |(prev, i, mut acc), curr| {
            acc.push((i, curr - prev));
            (curr + 1, i + 1, acc)
        }).2
    })
}

pub fn rep_combi(n: usize, r: usize) -> impl Iterator<Item = Vec<usize>> {
    rep_combi_count(n, r).map(|list| { 
        list.into_iter().fold(vec![], |mut acc, (i, counts)| { 
            acc.append(&mut vec![i; counts]);
            acc
        })
    })
}


#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn test_combi() { 
        let a = combi(5, 3);
        assert_eq!(a.count(), 10);
    }

    #[test]
    fn test_combi_zero() { 
        let a = combi(5, 0);
        assert_eq!(a.count(), 1);
    }

    #[test]
    fn test_combi_from_zero() { 
        let a = combi(0, 0);
        assert_eq!(a.count(), 1);

        let a = combi(0, 1);
        assert_eq!(a.count(), 0);
    }

    #[test]
    fn test_rep_combi() { 
        let a = rep_combi(5, 3);
        assert_eq!(a.count(), 35);
    }

    #[test]
    fn test_rep_combi_zero() { 
        let a = rep_combi(5, 0);
        assert_eq!(a.count(), 1);
    }

    #[test]
    fn test_rep_combi_from_zero() { 
        let a = rep_combi(0, 0);
        assert_eq!(a.count(), 1);

        let a = rep_combi(0, 1);
        assert_eq!(a.count(), 0);
    }
}