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
    nnz: Vec<Vec<Col>>, // [row -> [col]]
    row_head: Vec<Col>,     // [row -> col]
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

        PivotFinder{ shape, nnz, row_head, row_wght, col_wght, cands, pivots, max_pivot }
    }

    fn rows(&self) -> Row { 
        self.shape.0
    }

    fn cols(&self) -> Col { 
        self.shape.1
    }

    fn is_candidate(&self, i: Row, j: Col) -> bool { 
        self.cands[i].contains(&j)
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
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;
    use sprs::CsMat;
    use crate::{math::matrix::sparse::CsMatExt, utils::collections::hashset};
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
        assert_eq!(pf.row_head, vec![0,1,2,1,1,0]);
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
}