// Implementation based on:
// 
// "Parallel Sparse PLUQ Factorization modulo p", Charles Bouillaguet, Claire Delaplace, Marie-Emilie Voge.
// https://hal.inria.fr/hal-01646133/document
// 
// see also: SpaSM (Sparse direct Solver Modulo p)
// https://github.com/cbouilla/spasm

use std::slice::Iter;
use std::cmp::Ordering;
use std::collections::{HashSet, HashMap, VecDeque};
use itertools::Itertools;
use sprs::{CsMat, PermOwned};
use crate::math::traits::{Ring, RingOps};

pub const MAX_PIVOTS: usize = 300_000;

pub fn find_pivots<R>(a: &CsMat<R>) -> Vec<(usize, usize)>
where R: Ring, for<'x> &'x R: RingOps<R> {
    find_pivots_upto(a, MAX_PIVOTS)
}

pub fn find_pivots_upto<R>(a: &CsMat<R>, max_pivots: usize) -> Vec<(usize, usize)>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let mut pf = PivotFinder::new(a, max_pivots);

    pf.find_fl_pivots();
    pf.find_fl_col_pivots();
    pf.find_cycle_free_pivots();

    pf.result()
}

pub fn perms_by_pivots<R>(a: &CsMat<R>, pivs: &Vec<(usize, usize)>) -> (PermOwned, PermOwned)
where R: Ring, for<'x> &'x R: RingOps<R> {
    use std::collections::BTreeSet;
    fn perm(n: usize, v: Vec<usize>) -> PermOwned { 
        let mut set: BTreeSet<_> = (0..n).collect();
        let mut vec: Vec<usize> = vec![];
        for i in v { 
            vec.push(i);
            set.remove(&i);
        }
        for i in set { 
            vec.push(i);
        }
        let mut inv = vec![0; n];
        for (i, j) in vec.into_iter().enumerate() {
            inv[j] = i;
        }
        PermOwned::new(inv)
    }
    let (m, n) = a.shape();
    let (rows, cols) = pivs.iter().cloned().unzip();
    (perm(m, rows), perm(n, cols))
}

type Row = usize;
type Col = usize;
struct PivotFinder {
    shape: (usize, usize),
    entries: Vec<Vec<Col>>,     // [row -> [col]]
    row_wght: Vec<f32>,         // [row -> weight]
    col_wght: Vec<f32>,         // [col -> weight]
    cands: Vec<HashSet<Col>>,   // [row -> [col]]
    pivots: HashMap<Col, Row>,  // col -> row
    max_pivots: usize
}

#[allow(unused)]
impl PivotFinder { 
    fn new<R>(a: &CsMat<R>, max_pivots: usize) -> Self
    where R: Ring, for<'x> &'x R: RingOps<R> {
        let shape = a.shape();
        let (m, n) = shape;
        let mut entries = vec![vec![]; m];
        let mut row_head = vec![n; m];
        let mut row_wght = vec![0f32; m];
        let mut col_wght = vec![0f32; n];
        let mut cands = vec![HashSet::new(); m];
        let pivots = HashMap::new();

        for (r, (i, j)) in a.iter() { 
            if r.is_zero() { continue }
            entries[i].push(j);

            let w = 1f32; // TODO
            row_wght[i] += w;
            col_wght[j] += w;

            if row_head[i] == n || j < row_head[i] { 
                row_head[i] = j;
            }

            if r.is_unit() { 
                cands[i].insert(j);
            }
        }

        PivotFinder{ shape, entries, row_wght, col_wght, cands, pivots, max_pivots }
    }

    fn nrows(&self) -> Row { 
        self.shape.0
    }

    fn ncols(&self) -> Col { 
        self.shape.1
    }

    fn is_empty_row(&self, i: Row) -> bool { 
        self.entries[i].is_empty()
    }

    fn row_head(&self, i: Row) -> Option<Col> {
        self.entries[i].first().copied()
    }

    fn nz_cols(&self, i: Row) -> Iter<Col> {
        self.entries[i].iter()
    }

    fn cmp_rows(&self, i1: Row, i2: Row) -> Ordering {
        if let Some(o) = self.row_wght[i1].partial_cmp(&self.row_wght[i2]) { 
            o.then(Ord::cmp(&i1, &i2))
        } else { 
            Ordering::Equal
        }
    }

    fn cmp_cols(&self, j1: Col, j2: Col) -> Ordering {
        if let Some(o) = self.col_wght[j1].partial_cmp(&self.col_wght[j2]) { 
            o.then(Ord::cmp(&j1, &j2))
        } else {
            Ordering::Equal
        }
    }

    fn is_candidate(&self, i: Row, j: Col) -> bool { 
        self.cands[i].contains(&j)
    }

    fn npivs(&self) -> usize { 
        self.pivots.len()
    }

    fn is_piv_col(&self, j: Col) -> bool { 
        self.pivots.contains_key(&j)
    }

    fn piv_row(&self, j: Col) -> Option<Row>{ 
        self.pivots.get(&j).copied()
    }

    fn set_pivot(&mut self, i: Row, j: Col) {
        self.pivots.insert(j, i);
    }

    fn can_insert(&self) -> bool { 
        self.npivs() < self.max_pivots
    }

    fn remain_rows(&self) -> Vec<Row> { 
        let occ = self.pivots.values().collect::<HashSet<_>>();
        (0 .. self.nrows())
            .filter(|&i| 
                !occ.contains(&i) && !self.is_empty_row(i)
            )
            .sorted_by(|&i1, &i2| 
                self.cmp_rows(i1, i2)
            )
            .collect_vec()
    }

    fn occupied_cols(&self) -> HashSet<Col> {
        self.pivots.values().fold(HashSet::new(), |mut res, &i| {
            for &j in self.nz_cols(i) { 
                res.insert(j);
            }
            res
        })
    }

    fn find_fl_pivots(&mut self) {
        for i in self.remain_rows() {
            if !self.can_insert() { break }

            let Some(j) = self.row_head(i) else { continue };

            if !self.is_piv_col(j) && self.is_candidate(i, j) {
                self.set_pivot(i, j);
            }
        }
    }

    fn find_fl_col_pivots(&mut self) {
        let mut occ = self.occupied_cols();

        for i in self.remain_rows() { 
            if !self.can_insert() { break }

            let mut cands = vec![];

            for &j in self.nz_cols(i) { 
                if !occ.contains(&j) && self.is_candidate(i, j) {
                    cands.push(j);
                }
            }

            let Some(j) = cands.into_iter().sorted_by(|&j1, &j2| 
                self.cmp_cols(j1, j2)
            ).next() else { continue };

            self.set_pivot(i, j);

            for &j in self.nz_cols(i) { 
                occ.insert(j);
            }
        }
    }

    fn find_cycle_free_pivots(&mut self) {
        let n = self.ncols();
        let mut w = RowWorker::new(n);

        for i in self.remain_rows() { 
            if !self.can_insert() { break }

            if let Some(j) = self.cycle_free_pivot_in(i, &mut w) { 
                self.set_pivot(i, j);
            }
            w.clear();
        }
     }

    fn cycle_free_pivot_in(&self, i: usize, w: &mut RowWorker) -> Option<Col> {
        
        //       j          j2
        //  i [  o   #   #   #      # ]    *: pivot,
        //       |           :             #: candidate,
        //       V           : rmv         o: queued,
        // i2 [  * --> o --> x        ]    x: entry
        //             |            
        //             V            
        //    [        * ------> o    ]
        //                       |  

        let mut queue = VecDeque::new();

        for &j in self.nz_cols(i) {
            if self.is_piv_col(j) {
                queue.push_back(j);
                w.set_occupied(j);
            } else if self.is_candidate(i, j) {
                w.set_candidate(j);
            }
        }

        while !queue.is_empty() && w.has_candidate() { 
            let j = queue.pop_front().unwrap();
            let i2 = self.piv_row(j).unwrap();

            for &j2 in self.nz_cols(i2) { 
                if self.is_piv_col(j2) && !w.is_occupied(j2) { 
                    queue.push_back(j2);
                    w.set_occupied(j2);
                } else if w.is_candidate(j2) { 
                    w.set_occupied(j2);
                    if !w.has_candidate() { break }
                }
            }
        }

        w.collect_candidates().into_iter().sorted_by(|&j1, &j2| 
            self.cmp_cols(j1, j2)
        ).next()
    }

    fn result(&self) -> Vec<(usize, usize)> { 
        use topological_sort::TopologicalSort;
        let mut ts = TopologicalSort::new();
        
        for (&j, &i) in self.pivots.iter() { 
            ts.insert(j);
            for &j2 in self.nz_cols(i) { 
                if j != j2 && self.is_piv_col(j2) {
                    ts.add_dependency(j, j2);
                }
            }
        }

        ts.into_iter().map(|j| {
            let i = self.pivots[&j];
            (i, j)
        }).collect_vec()
    }
}

struct RowWorker { 
    status: Vec<i8>,
    ncand: usize
}

impl RowWorker {
    fn new(size: usize) -> Self { 
        let status = vec![0; size];
        RowWorker { status, ncand: 0 }
    }

    fn clear(&mut self) {
        self.status.fill(0);
        self.ncand = 0;
    }

    fn has_candidate(&self) -> bool { 
        self.ncand > 0
    }

    fn is_candidate(&self, i: usize) -> bool { 
        self.status[i] == 1
    }

    fn set_candidate(&mut self, i: usize) { 
        self.status[i] = 1;
        self.ncand += 1;
    }

    fn collect_candidates(&self) -> Vec<usize> { 
        self.status.iter().enumerate().filter_map(|(i, &s)| 
            if s == 1 { 
                Some(i)
            } else { 
                None
            }
        ).collect()
    }

    fn is_occupied(&self, i: usize) -> bool { 
        self.status[i] == -1
    }

    fn set_occupied(&mut self, i: usize) { 
        if self.is_candidate(i) { 
            self.ncand -= 1;
        }
        self.status[i] = -1;
    }
}

#[cfg(test)]
mod tests {
    use std::collections::{HashSet, HashMap};
    use num_traits::{Zero, One};
    use sprs::CsMat;
    use crate::math::matrix::{sparse::CsMatExt, pivot::perms_by_pivots};
    use crate::utils::collections::{hashmap, hashset};
    use super::*;
 
    #[test]
    fn init() {
        let a = CsMat::csc_from_vec((6, 9), vec![
            1, 0, 1, 0, 0, 1, 1, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 2, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 1, 0, 3, 0, 0, 0, 0,
            0, 1, 0, 1, 0, 0, 1, 0, 1,
            1, 0, 1, 0, 1, 1, 0, 1, 1
        ]);
        let pf = PivotFinder::new(&a, MAX_PIVOTS);

        assert_eq!(pf.entries, vec![
            vec![0,2,5,6,8], 
            vec![1,2,3,5,7], 
            vec![2,3,7,8], 
            vec![1,2,4], 
            vec![1,3,6,8], 
            vec![0,2,4,5,7,8]]
        );
        assert_eq!(pf.row_wght, vec![5f32,5f32,4f32,3f32,4f32,6f32]);
        assert_eq!(pf.col_wght, vec![2f32,3f32,5f32,3f32,2f32,3f32,2f32,3f32,4f32]);
        assert_eq!(pf.cands, vec![
            hashset!{0,2,5,6,8},
            hashset!{1,2,3,5},
            hashset!{2,3,7,8},
            hashset!{1,2},
            hashset!{1,3,6,8},
            hashset!{0,2,4,5,7,8}
        ]);
    }

    #[test]
    fn row_head() {
        let a = CsMat::csc_from_vec((4, 4), vec![
            1, 0, 1, 0,
            0, 1, 1, 1,
            0, 0, 0, 0,
            0, 0, 1, 1,
        ]);
        let pf = PivotFinder::new(&a, MAX_PIVOTS);

        assert_eq!(pf.row_head(0), Some(0));
        assert_eq!(pf.row_head(1), Some(1));
        assert_eq!(pf.row_head(2), None);
        assert_eq!(pf.row_head(3), Some(2));
    }

    #[test]
    fn rows_cols() {
        let a = CsMat::<i32>::csc_from_vec((4, 3), vec![]);
        let pf = PivotFinder::new(&a, MAX_PIVOTS);
        assert_eq!(pf.nrows(), 4);
        assert_eq!(pf.ncols(), 3);
    }

    #[test]
    fn remain_rows() {
        let a = CsMat::csc_from_vec((4, 4), vec![
            1, 0, 1, 0,
            0, 1, 1, 1,
            0, 0, 0, 0,
            0, 0, 1, 1,
        ]);
        let mut pf = PivotFinder::new(&a, MAX_PIVOTS);

        assert_eq!(pf.remain_rows(), vec![0,3,1]);

        pf.set_pivot(0, 0);

        assert_eq!(pf.remain_rows(), vec![3,1]);

        pf.set_pivot(1, 1);

        assert_eq!(pf.remain_rows(), vec![3]);
    }

    #[test]
    fn occupied_cols() {
        let a = CsMat::csc_from_vec((4, 4), vec![
            1, 0, 1, 0,
            0, 1, 1, 1,
            0, 0, 0, 0,
            0, 0, 1, 1,
        ]);
        let mut pf = PivotFinder::new(&a, MAX_PIVOTS);

        assert_eq!(pf.occupied_cols(), hashset!{});

        pf.set_pivot(0, 0);

        assert_eq!(pf.occupied_cols(), hashset!{0,2});

        pf.set_pivot(1, 1);

        assert_eq!(pf.occupied_cols(), hashset!{0,1,2,3});
    }

    #[test]
    fn find_fl_pivots() {
        let a = CsMat::csc_from_vec((6, 9), vec![
            1, 0, 1, 0, 0, 1, 1, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 1, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 1, 1
        ]);
        let mut pf = PivotFinder::new(&a, MAX_PIVOTS);

        pf.find_fl_pivots();

        assert_eq!(pf.pivots, hashmap!{2 => 4, 1 => 3, 5 => 5, 0 => 0} );
    }

    #[test]
    fn find_fl_col_pivots() { 
        let a = CsMat::csc_from_vec((6, 9), vec![
            1, 0, 0, 0, 0, 1, 0, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 1, 0, 1, 0
        ]);
        let mut pf = PivotFinder::new(&a, MAX_PIVOTS);

        pf.find_fl_col_pivots();

        assert_eq!(pf.pivots, hashmap!{2 => 4, 4 => 3, 0 => 0, 7 => 5, 3 => 2} );
    }

    #[test]
    fn find_fl_row_col_pivots() { 
        let a = CsMat::csc_from_vec((6, 9), vec![
            1, 0, 0, 0, 0, 1, 0, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 1, 0, 1, 0
        ]);
        let mut pf = PivotFinder::new(&a, MAX_PIVOTS);

        pf.find_fl_pivots();

        assert_eq!(pf.pivots, hashmap!{0 => 0, 1 => 3, 2 => 4} );

        pf.find_fl_col_pivots();

        assert_eq!(pf.pivots, hashmap!{0 => 0, 1 => 3, 2 => 4, 7 => 5, 3 => 2} );
    }

    #[test]
    fn find_cycle_free_pivots() {
        let a = CsMat::csc_from_vec((6, 9), vec![
            1, 0, 0, 0, 0, 1, 0, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 1, 0, 1, 0
        ]);
        let mut pf = PivotFinder::new(&a, MAX_PIVOTS);

        pf.find_cycle_free_pivots();

        assert_eq!(pf.pivots, hashmap!{2 => 4, 4 => 3, 0 => 0, 1 => 5, 3 => 2});
    }

    #[test]
    fn zero() { 
        let a = CsMat::csc_from_vec((1, 1), vec![
            0
        ]);
        let pivs = find_pivots(&a);
        let r = pivs.len();
        assert_eq!(r, 0);
    }

    #[test]
    fn id_1() { 
        let a = CsMat::csc_from_vec((1, 1), vec![
            1
        ]);
        let pivs = find_pivots(&a);
        let r = pivs.len();
        assert_eq!(r, 1);
    }

    #[test]
    fn id_2() { 
        let a = CsMat::csc_from_vec((2, 2), vec![
            1, 0, 0, 1
        ]);
        let pivs = find_pivots(&a);
        let r = pivs.len();
        assert_eq!(r, 2);
    }

    #[test]
    fn result() { 
        let a = CsMat::csc_from_vec((6, 9), vec![
            1, 0, 0, 0, 0, 1, 0, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 1, 0, 1, 0
        ]);
        let pivs = find_pivots(&a);
        let r = pivs.len();
        assert_eq!(r, 5);
        
        let (p, q) = perms_by_pivots(&a, &pivs);
        let b = a.permute(p.view(), q.view()).to_dense();

        assert!((0..r).all(|i| b[[i, i]].is_one()));
        assert!((0..r).all(|j| {
            (j+1..r).all(|i| b[[i, j]].is_zero())
        }));
    }

    #[test]
    fn rand() {
        let d = 0.1;
        let shape = (60, 80);
        let a: CsMat<i32> = CsMat::rand(shape, d);

        let pivs = find_pivots(&a);
        let r = pivs.len();
        assert!(r > 10);
        
        let (p, q) = perms_by_pivots(&a, &pivs);
        let b = a.permute(p.view(), q.view()).to_dense();

        assert!((0..r).all(|i| b[[i, i]].is_one()));
        assert!((0..r).all(|j| {
            (j+1..r).all(|i| b[[i, j]].is_zero())
        }))
    }
}