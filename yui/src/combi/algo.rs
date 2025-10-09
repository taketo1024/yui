use crate::combi::Partition;

// ref: https://en.wikipedia.org/wiki/Littlewood–Richardson_rule
pub fn lr_num(l: &Partition, m: &Partition, n: &Partition) -> usize { 
    if l.sum() + m.sum() != n.sum() { return 0 }
    if !n.contains(l) || !n.contains(m) { return 0 }

    fn is_target(l: &Partition, n: &Partition, i: usize, j: usize) -> bool {
        !l.contains_pt(i, j) && n.contains_pt(i, j)
    }

    fn helper(
        l: &Partition,
        m: &Partition,
        n: &Partition,
        i: usize,
        j: usize,
        tab: &mut Vec<Vec<usize>>,
        w: &mut Vec<usize>,
        count: &mut usize,
    ) {
        // reached final row
        if i == n.len() {
            if is_lattice_word(w) && cmp(m, w).is_eq() {
                *count += 1;
            }
            return;
        }

        // (i, j) not available in n/l
        if !is_target(l, n, i, j) { 
            // not yet reached to the right-end of l?
            if j > 0 && j >= l[i] { 
                // if so, move left
                helper(l, m, n, i, j - 1, tab, w, count);
            } else { 
                // otherwise, move down
                helper(l, m, n, i + 1, j, tab, w, count);
            }
            return;
        }

        // rows must be weakly increasing
        let max = if is_target(l, n, i, j + 1) { 
            tab[i][j + 1]
        } else { 
            // restriction from lattice-word condition.
            w.iter().max().map(|x| x + 1).unwrap_or(0)
        };

        // columns must be strictly increasing
        let min = if i > 0 && is_target(l, n, i - 1, j) { 
            tab[i - 1][j] + 1
        } else { 
            0
        };

        // try all possible numbers for (i, j)
        for x in min..=max { 
            w.push(x);

            if is_lattice_word(w) && cmp(m, w).is_ge() { 
                tab[i][j] = x;
                
                // to next depth
                if j > 0 && is_target(l, n, i, j - 1) { 
                    helper(l, m, n, i, j - 1, tab, w, count);
                } else { 
                    helper(l, m, n, i + 1, n[i + 1], tab, w, count);
                }
            }

            w.pop();
        }
    }

    let mut tab = vec![vec![0; n[0]]; n.len()]; // rectangle fitting n. 
    let mut w = Vec::with_capacity(m.sum());
    let mut count = 0;

    helper(l, m, n, 0, n[0], &mut tab, &mut w, &mut count);

    count
}

// ref: https://en.wikipedia.org/wiki/Lattice_word
pub fn is_lattice_word(w: &[usize]) -> bool { 
    let Some(m) = w.iter().max().copied() else { 
        return true
    };
    let mut counts = vec![0; m + 1];

    for &x in w.iter() {
        counts[x] += 1;
        if x > 0 && counts[x - 1] < counts[x] { 
            return false;
        }
    }
    true
}

// lex-order
fn cmp(m: &Partition, w: &Vec<usize>) -> std::cmp::Ordering {
    use std::cmp::Ordering::*;

    let Some(&max) = w.iter().max() else {
        return usize::cmp(&m.len(), &0)
    };

    for i in 0..=max { 
        let c = w.iter().filter(|&&j| i == j).count();
        if m[i] != c { 
            return usize::cmp(&m[i], &c)
        }
    }

    if m.len() > max + 1 { 
        Greater
    } else { 
        Equal
    }
}

#[cfg(test)]
mod tests {
    use num_traits::Zero;

    use crate::poly::PolyN;
    use crate::AddMon;

    use super::*;

    #[test]
    fn test_is_lattice_empty() {
        assert!(is_lattice_word(&[]));
    }

    #[test]
    fn test_is_lattice_single_element() {
        assert!(is_lattice_word(&[0]));
        assert!(!is_lattice_word(&[1]));
    }

    #[test]
    fn test_is_lattice_word() {
        assert!( is_lattice_word(&[0,0,0,1,1,0,1,0]));
        assert!(!is_lattice_word(&[0,1,0,1,1,0,0,0]));
    }

    #[test]
    fn test_cmp() { 
        use std::cmp::Ordering::*;

        assert_eq!(cmp(&Partition::from([]), &vec![]), Equal);
        assert_eq!(cmp(&Partition::from([1]), &vec![0]), Equal);
        assert_eq!(cmp(&Partition::from([2]), &vec![0, 0]), Equal);
        assert_eq!(cmp(&Partition::from([1,1]), &vec![0, 1]), Equal);
        assert_eq!(cmp(&Partition::from([2,1]), &vec![0, 0, 1]), Equal);
        assert_eq!(cmp(&Partition::from([2,1,1]), &vec![0, 0, 1]), Greater);
        assert_eq!(cmp(&Partition::from([2,1]), &vec![0, 0, 1, 2]), Less);
    }

    #[test]
    fn test_lr_num_zero_sum() {
        let l = Partition::from([]);
        let m = Partition::from([]);
        let n = Partition::from([]);
        assert_eq!(lr_num(&l, &m, &n), 1);
    }

    #[test]
    fn test_lr_num_sum_mismatch() {
        let l = Partition::from([1]);
        let m = Partition::from([1]);
        let n = Partition::from([1]);
        assert_eq!(lr_num(&l, &m, &n), 0);
    }

    #[test]
    fn test_lr_num_n_does_not_contain_l() {
        let l = Partition::from([2]);
        let m = Partition::from([1]);
        let n = Partition::from([1, 1, 1]);
        assert_eq!(lr_num(&l, &m, &n), 0);
    }

    #[test]
    fn test_lr_num_n_does_not_contain_m() {
        let l = Partition::from([1]);
        let m = Partition::from([2]);
        let n = Partition::from([1, 1, 1]);
        assert_eq!(lr_num(&l, &m, &n), 0);
    }

    #[test]
    fn test_lr_num_simple_case_1() {
        let l = Partition::from([1]);
        let m = Partition::from([1]);
        let n = Partition::from([2]);
        // |.|0|
        assert_eq!(lr_num(&l, &m, &n), 1); 
    }

    #[test]
    fn test_lr_num_simple_case_2() {
        let l = Partition::from([1, 0]);
        let m = Partition::from([1, 1]);
        let n = Partition::from([2, 1]);
        // |.|0|
        // |1|
        assert_eq!(lr_num(&l, &m, &n), 1);
    }

    #[test]
    fn test_lr_num_simple_case_3() {
        let l = Partition::from([1, 0]);
        let m = Partition::from([2]);
        let n = Partition::from([2, 1]);
        // |.|0|
        // |0|
        assert_eq!(lr_num(&l, &m, &n), 1);
    }

    #[test]
    fn test_lr_num_1() {
        let l = Partition::from([3,3,1]);
        let m = Partition::from([5,2]);
        let n = Partition::from([5,4,3,2]);
        assert_eq!(lr_num(&l, &m, &n), 1);
    }

    #[test]
    fn test_lr_num_2() {
        let l = Partition::from([3,3,1]);
        let m = Partition::from([4,3]);
        let n = Partition::from([5,4,3,2]);
        assert_eq!(lr_num(&l, &m, &n), 2);
    }

    #[test]
    fn test_lr_num_3() {
        let l = Partition::from([3,3,1]);
        let m = Partition::from([4,2,1]);
        let n = Partition::from([5,4,3,2]);
        assert_eq!(lr_num(&l, &m, &n), 3);
    }

    #[test]
    fn test_lr_num_4() {
        let l = Partition::from([3,3,1]);
        let m = Partition::from([4,1,1,1]);
        let n = Partition::from([5,4,3,2]);
        assert_eq!(lr_num(&l, &m, &n), 1);
    }

    #[test]
    fn test_lr_num_5() {
        let l = Partition::from([2,1]);
        let m = Partition::from([3,2]);
        let n = Partition::from([4,4]);
        assert_eq!(lr_num(&l, &m, &n), 1);
    }

    #[test]
    fn test_schur_poly_eqn_1() { 
        type P = PolyN<'x', i64>;

        let n = 3;
        let p = Partition::from([2,1]);
        let q = Partition::from([2,1]);
        let s_p = P::schur_poly(n, &p);
        let s_q = P::schur_poly(n, &q);
        let val = s_p * s_q;

        let d = p.sum() + q.sum();
        let expected = P::sum(
            Partition::all_partitions(d).map(|r|{
                let c = lr_num(&p, &q, &r) as i64;
                if c > 0 { 
                    P::from_const(c) * P::schur_poly(n, &r)
                } else { 
                    P::zero()
                }
            })
        );

        assert_eq!(val, expected)
    }
    
        #[test]
    fn test_schur_poly_eqn_2() { 
        type P = PolyN<'x', i64>;

        let n = 3;
        let p = Partition::from([2,1]);
        let q = Partition::from([3,2]);
        let s_p = P::schur_poly(n, &p);
        let s_q = P::schur_poly(n, &q);
        let val = s_p * s_q;

        let d = p.sum() + q.sum();
        let expected = P::sum(
            Partition::all_partitions(d).map(|r|{
                let c = lr_num(&p, &q, &r) as i64;

                if c > 0 { 
                    let c = P::from_const(c);
                    let s_r = P::schur_poly(n, &r);
                    c * s_r
                } else { 
                    P::zero()
                }
            })
        );

        assert_eq!(val, expected)
    }
}
