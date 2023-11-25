// Implementation based on:
// 
// "Parallel Sparse PLUQ Factorization modulo p", Charles Bouillaguet, Claire Delaplace, Marie-Emilie Voge.
// https://hal.inria.fr/hal-01646133/document
// 
// see also: SpaSM (Sparse direct Solver Modulo p)
// https://github.com/cbouilla/spasm

use std::cell::RefCell;
use std::slice::Iter;
use std::cmp::Ordering;
use std::collections::VecDeque;
use std::sync::{Mutex, RwLock};
use ahash::AHashSet;
use itertools::Itertools;
use log::info;
use rayon::prelude::{IntoParallelRefIterator, ParallelIterator};
use sprs::PermOwned;
use thread_local::ThreadLocal;
use yui::{Ring, RingOps};
use yui::algo::top_sort;
use super::*;

const LOG_THRESHOLD: usize = 10_000;

#[derive(PartialEq, Eq, Clone, Copy)]
pub enum PivotType { 
    Rows, Cols
}

pub fn find_pivots<R>(a: &SpMat<R>, piv_type: PivotType) -> Vec<(usize, usize)>
where R: Ring, for<'x> &'x R: RingOps<R> {
    if a.is_zero() { 
        return vec![];
    }
    
    let mut pf = PivotFinder::new(a, piv_type);
    pf.find_pivots();
    pf.result()
}

pub fn perms_by_pivots<R>(a: &SpMat<R>, pivs: &[(usize, usize)]) -> (PermOwned, PermOwned)
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

pub struct PivotFinder {
    str: MatrixStr,
    pivots: PivotData,
    piv_type: PivotType
}

impl PivotFinder { 
    pub fn new<R>(a: &SpMat<R>, piv_type: PivotType) -> Self
    where R: Ring, for<'x> &'x R: RingOps<R> {
        let str = MatrixStr::new(a, piv_type);
        let pivots = PivotData::new(a, piv_type);
        PivotFinder{ str, pivots, piv_type }
    }

    pub fn find_pivots(&mut self) {
        if self.should_report() {
            info!("find pivots: {:?}", self.str.shape());
        }

        self.find_fl_pivots();
        self.find_fl_col_pivots();
        self.find_cycle_free_pivots();

        if self.should_report() {
            info!("found {} pivots.", self.pivots.count());
        }
    }

    pub fn result(&self) -> Vec<(usize, usize)> { 
        let is_row_type = self.piv_type == PivotType::Rows;

        let tree = self.pivots.iter().map(|(i, j)| { 
            let list = self.str.cols_in(i).filter(|&&j2|
                j != j2 && self.pivots.has_col(j2)
            ).copied().collect_vec();
            (j, list)
        });
        
        let sorted = top_sort(tree).unwrap();

        sorted.into_iter().map(|j| {
            let i = self.pivots.row_for(j).unwrap();
            if is_row_type { (i, j) } else { (j, i) }
        }).collect_vec()
    }

    fn rows(&self) -> Row { 
        self.str.shape.0
    }

    fn cols(&self) -> Col { 
        self.str.shape.1
    }

    fn remain_rows(&self) -> impl Iterator<Item = Row> { 
        let piv_rows: AHashSet<_> = self.pivots.iter().map(|(i, _)| i).collect();
        let m = self.rows();

        (0 .. m).filter(|&i| 
            !piv_rows.contains(&i) && !self.str.is_empty_row(i)
        ).sorted_by(|&i1, &i2| 
            self.str.cmp_rows(i1, i2)
        )
    }

    fn occupied_cols(&self) -> AHashSet<Col> {
        self.pivots.iter().fold(AHashSet::new(), |mut res, (i, _)| {
            for &j in self.str.cols_in(i) { 
                res.insert(j);
            }
            res
        })
    }

    fn find_fl_pivots(&mut self) {
        let remain_rows: Vec<_> = self.remain_rows().collect();

        for i in remain_rows {
            let Some(j) = self.str.head_col_in(i) else { continue };

            if !self.pivots.has_col(j) && self.str.is_candidate(i, j) {
                self.pivots.set(i, j);
            }
        }

        let piv_count = self.pivots.count();

        if self.should_report() { 
            info!("  fl-pivots: +{}.", piv_count);
        }
    }

    fn find_fl_col_pivots(&mut self) {
        let before_piv_count = self.pivots.count();

        let remain_rows: Vec<_> = self.remain_rows().collect();
        let mut occ_cols = self.occupied_cols();

        for i in remain_rows { 
            let mut cands = vec![];

            for &j in self.str.cols_in(i) { 
                if !occ_cols.contains(&j) && self.str.is_candidate(i, j) {
                    cands.push(j);
                }
            }

            let Some(j) = cands.into_iter().sorted_by(|&j1, &j2| 
                self.str.cmp_cols(j1, j2)
            ).next() else { continue };

            self.pivots.set(i, j);

            for &j in self.str.cols_in(i) { 
                occ_cols.insert(j);
            }
        }

        let piv_count = self.pivots.count();
        
        if self.should_report() { 
            info!("  fl-col-pivots: +{}, total: {}.", piv_count - before_piv_count, piv_count);
        }
    }

    fn find_cycle_free_pivots(&mut self) {
        let before_piv_count = self.pivots.count();

        if crate::config::is_multithread_enabled() { 
            self.find_cycle_free_pivots_m();
        } else {
            self.find_cycle_free_pivots_s();
        }

        let piv_count = self.pivots.count();

        if self.should_report() { 
            info!("  cycle-free-pivots: +{}, total: {}.", piv_count - before_piv_count, piv_count);
        }
    }

    fn find_cycle_free_pivots_s(&mut self) {
        let remain_rows: Vec<_> = self.remain_rows().collect();
        let total_rows = remain_rows.len();

        let mut row_count = 0;
        let mut piv_count = self.pivots.count();

        if self.should_report() { 
            info!("  start find-cycle-free-pivots: {total_rows} rows");
        }

        let n = self.cols();
        let mut w = RowWorker::new(n);

        for i in remain_rows { 
            if let Some(j) = w.find_cycle_free_pivots(i, &self.str, &self.pivots) { 
                self.pivots.set(i, j);
            }

            if self.should_report() {
                row_count += 1;
                if row_count % LOG_THRESHOLD == 0 { 
                    let c = self.pivots.count();
                    info!("    [{row_count}/{total_rows}] +{}, total: {}.", c - piv_count, c);
                    piv_count = c;
                }
            }
        }
     }

     fn find_cycle_free_pivots_m(&mut self) {
        let n = self.cols();
        let remain_rows = self.remain_rows().collect_vec();
        let total_rows = remain_rows.len();

        let pivots = RwLock::new(
            std::mem::take(&mut self.pivots)
        );
        let loc_pivots_tls = ThreadLocal::new();
        let loc_worker_tls = ThreadLocal::new();

        let (row_count, piv_count) = if self.should_report() { 
            let n_threads = std::thread::available_parallelism().map(|x| x.get()).unwrap_or(1);

            info!("  start find-cycle-free-pivots: {total_rows} rows (multi-thread: {n_threads})");

            let row_count = Mutex::new(0);
            let piv_count = Mutex::new(self.pivots.count());

            (Some(row_count), Some(piv_count))
        } else { 
            (None, None)
        };

        remain_rows.par_iter().for_each(|&i| { 
            let mut loc_pivots = init_tls(&loc_pivots_tls, || 
                pivots.read().unwrap().clone()
            ).borrow_mut();

            let mut w = init_tls(&loc_worker_tls, || 
                RowWorker::new(n)
            ).borrow_mut();

            loc_pivots.update_from(&pivots.read().unwrap());
            w.init(i, &self.str, &loc_pivots);

            self.find_cycle_free_pivots_in(&pivots, &mut loc_pivots, &mut w);

            if self.should_report() { 
                let row_count = {
                    let mut row_count = row_count.as_ref().unwrap().lock().unwrap();
                    *row_count += 1;
                    *row_count
                };
            
                if row_count % 10_000 == 0 { 
                    let mut piv_count = piv_count.as_ref().unwrap().lock().unwrap();
                    let c = loc_pivots.count();
                    info!("    [{row_count}/{total_rows}] +{}, total: {}.", c - *piv_count, c);
                    *piv_count = c;
                }
            }
        });

        self.pivots = pivots.into_inner().unwrap();
     }

     #[inline(never)] // for profilability
     fn find_cycle_free_pivots_in(&self, pivots: &RwLock<PivotData>, loc_pivots: &mut PivotData, w: &mut RowWorker) {
        loop { 
            w.traverse(&self.str, loc_pivots);
    
            let Some(j) = w.choose_candidate(&self.str) else {
                break
            };
            
            // If changes are made in other threads, update `loc_pivots` and retry.
            // Otherwise, modify `pivots` and exit.
        
            let mut pivots = pivots.write().unwrap();
            w.update_diff(&loc_pivots, &pivots);
            
            if w.should_retry() { 
                loc_pivots.update_from(&pivots);
                continue
            } else { 
                pivots.set(w.row, j);
                break
            }    
        }
     }

     fn should_report(&self) -> bool { 
        self.rows() > LOG_THRESHOLD && log::max_level() >= log::LevelFilter::Info
     }
}

fn init_tls<T, F>(tl: &ThreadLocal<RefCell<T>>, f: F) -> &RefCell<T>
where T: Send, F: FnOnce() -> T {
    tl.get_or(|| RefCell::new( f() ) )
}

struct MatrixStr { 
    shape: (usize, usize),
    entries: Vec<Vec<Col>>,     // [row -> [col]]
    cands: Vec<AHashSet<Col>>,  // [row -> [col]]
    row_wght: Vec<f32>,         // [row -> weight]
    col_wght: Vec<f32>,         // [col -> weight]
}

impl MatrixStr { 
    fn new<R>(a: &SpMat<R>, piv_type: PivotType) -> Self
    where R: Ring, for<'x> &'x R: RingOps<R> { 
        let a_view = a.view();
        let a = if piv_type == PivotType::Rows { 
            a_view 
        } else { 
            a_view.transpose()
        };

        let shape = a.shape();
        let (m, n) = shape;
        let mut entries = vec![vec![]; m];
        let mut row_wght = vec![0f32; m];
        let mut col_wght = vec![0f32; n];
        let mut cands = vec![AHashSet::new(); m];

        for (i, j, r) in a.iter() { 
            assert!(!r.is_zero());
            
            entries[i].push(j);

            let w = 1f32; // TODO
            row_wght[i] += w;
            col_wght[j] += w;

            if r.is_unit() { 
                cands[i].insert(j);
            }
        }

        Self { shape, entries, cands, row_wght, col_wght }
    }

    fn shape(&self) -> (usize, usize) {
        self.shape
    }

    fn is_empty_row(&self, i: Row) -> bool { 
        self.entries[i].is_empty()
    }

    fn head_col_in(&self, i: Row) -> Option<Col> {
        self.entries[i].first().copied()
    }

    fn cols_in(&self, i: Row) -> Iter<Col> {
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
}

#[derive(Clone, Default)]
struct PivotData { 
    data: Vec<Option<Row>>,   // col -> row
    indices: Vec<Col>
}

impl PivotData { 
    fn new<R>(a: &SpMat<R>, piv_type: PivotType) -> Self
    where R: Ring, for<'x> &'x R: RingOps<R> { 
        let n = if piv_type == PivotType::Rows { 
            a.cols()
        } else { 
            a.rows()
        };
        let data = vec![None; n];
        let indices = vec![];
        Self { data, indices }
    }

    fn count(&self) -> usize { 
        self.indices.len()
    }

    fn has_col(&self, j: Col) -> bool { 
        self.data[j].is_some()
    }

    fn row_for(&self, j: Col) -> Option<Row> { 
        self.data[j]
    }

    fn set(&mut self, i: Row, j: Col) {
        assert!(!self.has_col(j));
        self.data[j] = Some(i);
        self.indices.push(j);
    }

    fn iter(&self) -> impl Iterator<Item = (Row, Col)> + '_ { 
        self.indices.iter().map(|&j| {
            let i = self.data[j].unwrap();
            (i, j)
        })
    }

    fn pivot_at(&self, k: usize) -> (Row, Col) { 
        let j = self.indices[k];
        let i = self.data[j].unwrap();
        (i, j)
    }

    fn update_from(&mut self, from: &Self) { 
        debug_assert!(self.count() <= from.count());
        for k in self.count() .. from.count() { 
            let (i, j) = from.pivot_at(k);
            self.set(i, j);
        }
    }
}

struct RowWorker { 
    row: usize,
    status: Vec<i8>,
    ncand: usize,
    queue: VecDeque<Col>,
    queued: AHashSet<Col>
}

impl RowWorker {
    fn new(size: usize) -> Self { 
        let status = vec![0; size];
        let queue = VecDeque::new();
        let queued = AHashSet::new();
        RowWorker {row: 0, status, ncand: 0, queue, queued }
    }

    fn clear(&mut self) {
        self.row = 0;
        self.status.fill(0);
        self.ncand = 0;
        self.queue.clear();
        self.queued.clear();
    }

    //  i [  o       #     # ]     [  o   x   x      # ]     [  o   x   x   x  # ]
    //    [  |               ]     [  |   :   :        ]     [  |   :   :   :    ]
    //    [  *   .   .       ] ~~> [  * - o - .        ] ~~> [  * - o - .   :    ]
    //    [                  ]     [      |            ]     [      |       :    ]
    //    [      *       .   ]     [      *       .    ]     [      *-------.    ]
    //
    //  o: queued, #: candidate, x: occupied

    fn find_cycle_free_pivots(&mut self, i: usize, str: &MatrixStr, pivots: &PivotData) -> Option<Col> {
        self.init(i, str, pivots);
        self.traverse(str, pivots);
        self.choose_candidate(str)
    }

    fn init(&mut self, i: usize, str: &MatrixStr, pivots: &PivotData) { 
        self.clear();
        self.row = i;

        for &j in str.cols_in(i) {
            if pivots.has_col(j) {
                self.enqueue(j);
                self.set_occupied(j);
            } else if str.is_candidate(i, j) {
                self.set_candidate(j);
            } else { 
                self.set_occupied(j);
            }
        }
    }

    fn traverse(&mut self, str: &MatrixStr, pivots: &PivotData) {
        if !self.has_candidate() { 
            return
        }

        while let Some(j) = self.dequeue() { 
            let i2 = pivots.row_for(j).unwrap();

            for &j2 in str.cols_in(i2) { 
                if pivots.has_col(j2) && !self.is_queued(j2) { 
                    self.enqueue(j2);
                }

                self.set_occupied(j2);

                if !self.has_candidate() { 
                    break 
                }
            }
        }
    }

    fn choose_candidate(&self, str: &MatrixStr) -> Option<Col> { 
        let n = self.status.len();
        (0 .. n)
            .filter(|&j| self.is_candidate(j))
            .sorted_by(|&j1, &j2| 
                str.cmp_cols(j1, j2)
            ).next()
        }

    fn update_diff(&mut self, loc_pivots: &PivotData, pivots: &PivotData) {
        debug_assert!(loc_pivots.count() <= pivots.count());
        for k in loc_pivots.count()..pivots.count() { 
            let j = pivots.indices[k];
            if self.is_candidate(j) || self.is_occupied(j) {
                self.enqueue(j);
                self.set_occupied(j);
            }
        }
    }

    fn should_retry(&self) -> bool { 
        !self.queue.is_empty()
    }

    fn has_candidate(&self) -> bool { 
        self.ncand > 0
    }

    fn is_candidate(&self, i: usize) -> bool { 
        self.status[i] == 1
    }

    fn set_candidate(&mut self, i: usize) { 
        assert_eq!(self.status[i], 0);
        self.status[i] = 1;
        self.ncand += 1;
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

    fn enqueue(&mut self, i: Col) { 
        self.queue.push_back(i);
        self.queued.insert(i);
    }

    fn dequeue(&mut self) -> Option<Col> { 
        self.queue.pop_front()
    }

    fn is_queued(&self, i: Col) -> bool { 
        self.queued.contains(&i)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::{Zero, One};
 
    #[test]
    fn str_init() {
        let a = SpMat::from_vec((6, 9), vec![
            1, 0, 1, 0, 0, 1, 1, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 2, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 1, 0, 3, 0, 0, 0, 0,
            0, 1, 0, 1, 0, 0, 1, 0, 1,
            1, 0, 1, 0, 1, 1, 0, 1, 1
        ]);
        let str = MatrixStr::new(&a, PivotType::Rows);

        assert_eq!(str.entries, vec![
            vec![0,2,5,6,8], 
            vec![1,2,3,5,7], 
            vec![2,3,7,8], 
            vec![1,2,4], 
            vec![1,3,6,8], 
            vec![0,2,4,5,7,8]]
        );
        assert_eq!(str.row_wght, vec![5f32,5f32,4f32,3f32,4f32,6f32]);
        assert_eq!(str.col_wght, vec![2f32,3f32,5f32,3f32,2f32,3f32,2f32,3f32,4f32]);
        assert_eq!(str.cands, vec![
            AHashSet::from_iter([0,2,5,6,8]),
            AHashSet::from_iter([1,2,3,5]),
            AHashSet::from_iter([2,3,7,8]),
            AHashSet::from_iter([1,2]),
            AHashSet::from_iter([1,3,6,8]),
            AHashSet::from_iter([0,2,4,5,7,8])
        ]);
    }

    #[test]
    fn str_row_head() {
        let a = SpMat::from_vec((4, 4), vec![
            1, 0, 1, 0,
            0, 1, 1, 1,
            0, 0, 0, 0,
            0, 0, 1, 1,
        ]);
        let str = MatrixStr::new(&a, PivotType::Rows);

        assert_eq!(str.head_col_in(0), Some(0));
        assert_eq!(str.head_col_in(1), Some(1));
        assert_eq!(str.head_col_in(2), None);
        assert_eq!(str.head_col_in(3), Some(2));
    }

    #[test]
    fn rows_cols() {
        let a = SpMat::<i32>::from_vec((4, 3), vec![]);
        let pf = PivotFinder::new(&a, PivotType::Rows);
        assert_eq!(pf.rows(), 4);
        assert_eq!(pf.cols(), 3);
    }

    #[test]
    fn pivot_data() { 
        let a = SpMat::from_vec((2, 4), vec![
            1, 0, 1, 0,
            0, 0, 1, 1,
        ]);
        let mut piv = PivotData::new(&a, PivotType::Rows);

        assert_eq!(piv.count(), 0);
        assert!(!piv.has_col(0));
        assert_eq!(piv.row_for(0), None);

        piv.set(1, 2);

        assert_eq!(piv.count(), 1);
        assert!(piv.has_col(2));
        assert_eq!(piv.row_for(2), Some(1));
    }

    #[test]
    fn remain_rows() {
        let a = SpMat::from_vec((4, 4), vec![
            1, 0, 1, 0,
            0, 1, 1, 1,
            0, 0, 0, 0,
            0, 0, 1, 1,
        ]);
        let mut pf = PivotFinder::new(&a, PivotType::Rows);

        assert_eq!(pf.remain_rows().collect_vec(), vec![0,3,1]);

        pf.pivots.set(0, 0);

        assert_eq!(pf.remain_rows().collect_vec(), vec![3,1]);

        pf.pivots.set(1, 1);

        assert_eq!(pf.remain_rows().collect_vec(), vec![3]);
    }

    #[test]
    fn pivots() {
        let a = SpMat::from_vec((4, 4), vec![
            1, 0, 1, 0,
            0, 1, 1, 1,
            0, 0, 0, 0,
            0, 0, 1, 1,
        ]);
        let mut pf = PivotFinder::new(&a, PivotType::Rows);

        assert_eq!(pf.occupied_cols(), AHashSet::new());

        pf.pivots.set(0, 0);

        assert_eq!(pf.occupied_cols(), AHashSet::from_iter([0,2]));

        pf.pivots.set(1, 1);

        assert_eq!(pf.occupied_cols(), AHashSet::from_iter([0,1,2,3]));
    }

    #[test]
    fn find_fl_pivots() {
        let a = SpMat::from_vec((6, 9), vec![
            1, 0, 1, 0, 0, 1, 1, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 1, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 1, 1
        ]);
        let mut pf = PivotFinder::new(&a, PivotType::Rows);

        pf.find_fl_pivots();

        assert_eq!(pf.pivots.iter().collect_vec(), vec![(4, 2), (3, 1), (5, 5), (0, 0)]);
    }

    #[test]
    fn find_fl_col_pivots() { 
        let a = SpMat::from_vec((6, 9), vec![
            1, 0, 0, 0, 0, 1, 0, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 1, 0, 1, 0
        ]);
        let mut pf = PivotFinder::new(&a, PivotType::Rows);

        pf.find_fl_col_pivots();

        assert_eq!(pf.pivots.iter().collect_vec(), vec![(4, 2), (3, 4), (0, 0), (5, 7), (2, 3)]);
    }

    #[test]
    fn find_fl_row_col_pivots() { 
        let a = SpMat::from_vec((6, 9), vec![
            1, 0, 0, 0, 0, 1, 0, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 1, 0, 1, 0
        ]);
        let mut pf = PivotFinder::new(&a, PivotType::Rows);

        pf.find_fl_pivots();

        assert_eq!(pf.pivots.iter().collect_vec(), vec![(4, 2), (3, 1), (0, 0)]);

        pf.find_fl_col_pivots();

        assert_eq!(pf.pivots.iter().collect_vec(), vec![(4, 2), (3, 1), (0, 0), (5, 7), (2, 3)]);
    }

    #[test]
    fn find_cycle_free_pivots_s() {
        let a = SpMat::from_vec((6, 9), vec![
            1, 0, 0, 0, 0, 1, 0, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 1, 0, 1, 0
        ]);
        let mut pf = PivotFinder::new(&a, PivotType::Rows);

        pf.find_cycle_free_pivots_s();

        assert_eq!(pf.pivots.iter().collect_vec(), vec![(4, 2), (3, 4), (0, 0), (5, 1), (2, 3)]);
    }

    #[test]
    fn find_cycle_free_pivots_m() {
        let a = SpMat::from_vec((6, 9), vec![
            1, 0, 0, 0, 0, 1, 0, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 1, 0, 1, 0
        ]);
        let mut pf = PivotFinder::new(&a, PivotType::Rows);

        pf.find_cycle_free_pivots_m();

        assert!(pf.pivots.count() >= 5);
    }

    #[test]
    fn zero() { 
        let a = SpMat::from_vec((1, 1), vec![
            0
        ]);
        let pivs = find_pivots(&a, PivotType::Rows);
        let r = pivs.len();
        assert_eq!(r, 0);
    }

    #[test]
    fn id_1() { 
        let a = SpMat::from_vec((1, 1), vec![
            1
        ]);
        let pivs = find_pivots(&a, PivotType::Rows);
        let r = pivs.len();
        assert_eq!(r, 1);
    }

    #[test]
    fn id_2() { 
        let a = SpMat::from_vec((2, 2), vec![
            1, 0, 0, 1
        ]);
        let pivs = find_pivots(&a, PivotType::Rows);
        let r = pivs.len();
        assert_eq!(r, 2);
    }

    #[test]
    fn result() { 
        let a = SpMat::from_vec((6, 9), vec![
            1, 0, 0, 0, 0, 1, 0, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 1, 0, 1, 0
        ]);
        let pivs = find_pivots(&a, PivotType::Rows);
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
    fn result_cols() { 
        let a = SpMat::from_vec((6, 9), vec![
            1, 0, 0, 0, 0, 1, 0, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 1, 0, 1, 0
        ]);
        let pivs = find_pivots(&a, PivotType::Cols);
        let r = pivs.len();
        assert_eq!(r, 6);
        
        let (p, q) = perms_by_pivots(&a, &pivs);
        let b = a.permute(p.view(), q.view()).to_dense();

        assert!((0..r).all(|i| b[[i, i]].is_one()));
        assert!((0..r).all(|i| {
            (i+1..r).all(|j| b[[i, j]].is_zero())
        }));
    }

    #[test]
    fn rand() {
        let d = 0.1;
        let shape = (60, 80);
        let a = SpMat::<i32>::rand(shape, d);

        let pivs = find_pivots(&a, PivotType::Rows);
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