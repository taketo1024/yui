use std::collections::{HashSet, HashMap};

use itertools::Itertools;
use sprs::CsMat;

use crate::math::traits::{Ring, RingOps};

pub fn find_pivots<R>(a: &CsMat<R>) -> Vec<(usize, usize)>
where R: Ring, for<'x> &'x R: RingOps<R> {
    let pf = PivotFinder::new(a);
    todo!()
}

type Row = usize;
type Col = usize;
struct PivotFinder {
    shape: (usize, usize),
    nnz: Vec<Vec<Col>>,     // [row -> [col]]
    row_wght: Vec<f32>,     // [row -> weight]
    col_wght: Vec<f32>,     // [col -> weight]
    cands: Vec<HashSet<Col>>, // [row -> [col]]
    pivots: HashMap<Col, Row>,  // col -> row
    max_pivot: usize
}

impl PivotFinder { 
    fn new<R>(a: &CsMat<R>) -> Self
    where R: Ring, for<'x> &'x R: RingOps<R> {
        let shape = a.shape();
        let (m, n) = shape;
        let mut nnz = vec![vec![]; m];
        let mut row_head = vec![n; m];
        let mut row_wght = vec![0f32; m];
        let mut col_wght = vec![0f32; n];
        let mut cands = vec![HashSet::new(); m];
        let pivots = HashMap::new();
        let max_pivot = 300000;

        for (r, (i, j)) in a.iter() { 
            if r.is_zero() { continue }
            nnz[i].push(j);

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

        PivotFinder{ shape, nnz, row_wght, col_wght, cands, pivots, max_pivot }
    }

    fn rows(&self) -> Row { 
        self.shape.0
    }

    fn cols(&self) -> Col { 
        self.shape.1
    }

    fn row_head(&self, i: Row) -> Option<Col> {
        self.nnz[i].first().map(|&j| j)
    }

    fn is_candidate(&self, i: Row, j: Col) -> bool { 
        self.cands[i].contains(&j)
    }

    fn npivs(&self) -> usize { 
        self.pivots.len()
    }

    fn has_pivot(&self, j: Col) -> bool { 
        self.pivots.contains_key(&j)
    }

    fn piv_row(&self, j: Col) -> Option<Row>{ 
        self.pivots.get(&j).copied()
    }

    fn set_pivot(&mut self, i: Row, j: Col) {
        self.pivots.insert(j, i);
    }

    fn remain_rows(&self) -> Vec<Row> { 
        let occ = self.pivots.values().collect::<HashSet<_>>();
        (0 .. self.rows()).into_iter()
            .filter(|i| !occ.contains(i))
            .sorted_by(|&i1, &i2| self.row_wght[i1].partial_cmp(&self.row_wght[i2]).unwrap())
            .collect_vec()
    }

    fn occupied_cols(&self) -> HashSet<Col> {
        self.pivots.values().fold(HashSet::new(), |mut res, &i| {
            for j in &self.nnz[i] { 
                res.insert(*j);
            }
            res
        })
    }

    fn find_fl_pivots(&mut self) {
        let m = self.rows();

        for i in 0..m {
            let Some(j) = self.row_head(i) else { continue };
            if !self.is_candidate(i, j) { continue }

            let should_ins = 
                self.piv_row(j).map(|i2| // occupied
                    self.row_wght[i] < self.row_wght[i2]
                ).unwrap_or(             // free
                    self.npivs() < self.max_pivot
                );

            if should_ins {
                self.set_pivot(i, j);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use std::collections::{HashSet, HashMap};
    use sprs::CsMat;
    use crate::{math::matrix::sparse::CsMatExt, utils::collections::*};
    use super::PivotFinder;
 
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
        let pf = PivotFinder::new(&a);

        assert_eq!(pf.nnz, vec![
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
        let pf = PivotFinder::new(&a);

        assert_eq!(pf.row_head(0), Some(0));
        assert_eq!(pf.row_head(1), Some(1));
        assert_eq!(pf.row_head(2), None);
        assert_eq!(pf.row_head(3), Some(2));
    }

    #[test]
    fn rows_cols() {
        let a = CsMat::<i32>::csc_from_vec((4, 3), vec![]);
        let pf = PivotFinder::new(&a);
        assert_eq!(a.rows(), 4);
        assert_eq!(a.cols(), 3);
    }

    #[test]
    fn remain_rows() {
        let a = CsMat::<i32>::csc_from_vec((4, 3), vec![]);
        let mut pf = PivotFinder::new(&a);

        assert_eq!(pf.remain_rows(), vec![0,1,2,3]);

        pf.set_pivot(0, 0);
        pf.set_pivot(2, 1);

        assert_eq!(pf.remain_rows(), vec![1,3]);
    }

    #[test]
    fn occupied_cols() {
        let a = CsMat::csc_from_vec((4, 4), vec![
            1, 0, 1, 0,
            0, 1, 1, 1,
            0, 0, 0, 0,
            0, 0, 1, 1,
        ]);
        let mut pf = PivotFinder::new(&a);

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
        let mut pf = PivotFinder::new(&a);

        pf.find_fl_pivots();

        assert_eq!(pf.pivots, hashmap!{0 => 0, 1 => 3, 2 => 4, 5 => 5} );
    }
}