// Implementation based on:
// 
// "Parallel Sparse PLUQ Factorization modulo p", Charles Bouillaguet, Claire Delaplace, Marie-Emilie Voge.
// https://hal.inria.fr/hal-01646133/document
// 
// see also: SpaSM (Sparse direct Solver Modulo p)
// https://github.com/cbouilla/spasm

use std::slice::Iter;
use std::cmp::Ordering;
use std::collections::VecDeque;
use ahash::AHashSet;
use itertools::Itertools;
use log::trace;
use sprs::PermOwned;

use yui::{Ring, RingOps};
use yui::algo::top_sort;
use super::*;
use super::util::perm_for_indices;

cfg_if::cfg_if! {
    if #[cfg(feature = "multithread")] {
        use std::cell::RefCell;
        use std::sync::RwLock;
        use thread_local::ThreadLocal;
        use rayon::prelude::*;
    }
}

const LOG_THRESHOLD: usize = 10_000;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum PivotType { 
    Rows, Cols
}

#[derive(Clone, Copy, Debug)]
pub enum PivotCondition { 
    One, Weight(f64), AnyUnit
}

impl PivotCondition { 
    fn is_cand<R>(&self, r: &R) -> bool 
    where R: Ring, for<'x> &'x R: RingOps<R> {
        match self {
            PivotCondition::One       => r.is_pm_one(),
            PivotCondition::Weight(w) => r.is_unit() && r.c_weight() <= *w,
            PivotCondition::AnyUnit   => r.is_unit(),
        }
    }
}

pub fn find_pivots<R>(a: &SpMat<R>, piv_type: PivotType, pivot_cond: PivotCondition) -> Vec<(usize, usize)>
where R: Ring, for<'x> &'x R: RingOps<R> {
    if a.is_zero() { 
        return vec![];
    }
    
    let mut pf = PivotFinder::new(a, piv_type, pivot_cond);
    pf.find_pivots();
    pf.result()
}

pub fn perms_by_pivots<R>(a: &SpMat<R>, pivs: &[(usize, usize)]) -> (PermOwned, PermOwned)
where R: Ring, for<'x> &'x R: RingOps<R> {
    let (m, n) = a.shape();
    (
        perm_for_indices(m, pivs.iter().map(|(i, _)| i)), 
        perm_for_indices(n, pivs.iter().map(|(_, j)| j))
    )
}

type Row = usize;
type Col = usize;

pub struct PivotFinder {
    str: MatrixStr,
    pivots: PivotData,
    piv_type: PivotType
}

impl PivotFinder { 
    pub fn new<R>(a: &SpMat<R>, piv_type: PivotType, pivot_cond: PivotCondition) -> Self
    where R: Ring, for<'x> &'x R: RingOps<R> {
        let str = MatrixStr::new(a, piv_type, pivot_cond);
        let pivots = PivotData::new(a, piv_type);
        PivotFinder{ str, pivots, piv_type }
    }

    pub fn find_pivots(&mut self) {
        trace!("pivots: {:?} ..", self.str.shape());

        self.find_fl_pivots();
        self.find_fl_col_pivots();
        self.find_cycle_free_pivots();

        trace!("pivots: {:?} => {}.", self.str.shape(), self.pivots.count());
    }

    pub fn result(&self) -> Vec<(usize, usize)> { 
        let tree = self.pivots.iter().map(|(i, j)| { 
            let list = self.str.cols_in(i).filter(|&&j2|
                j != j2 && self.pivots.has_col(j2)
            ).copied().collect_vec();
            (j, list)
        });
        
        let sorted = top_sort(tree).unwrap();
        let is_row_type = self.piv_type == PivotType::Rows;
        
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

        trace!("  fl-pivots: +{}.", piv_count);
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
        
        trace!("  fl-col-pivots: +{}, total: {}.", piv_count - before_piv_count, piv_count);
    }

    fn find_cycle_free_pivots(&mut self) {
        let before_piv_count = self.pivots.count();

        cfg_if::cfg_if! { 
            if #[cfg(feature = "multithread")] { 
                self.find_cycle_free_pivots_m();
            } else { 
                self.find_cycle_free_pivots_s();
            }
        }

        let piv_count = self.pivots.count();

        trace!("  cycle-free-pivots: +{}, total: {}.", piv_count - before_piv_count, piv_count);
    }

    #[allow(unused)]
    fn find_cycle_free_pivots_s(&mut self) {
        let remain_rows: Vec<_> = self.remain_rows().collect();
        let total_rows = remain_rows.len();

        trace!("  start find-cycle-free-pivots: {total_rows} rows");

        let n = self.cols();
        let mut w = RowWorker::new(n);
        let mut row_count = 0;

        for i in remain_rows { 
            if let Some(j) = w.find_cycle_free_pivots(i, &self.str, &self.pivots) { 
                self.pivots.set(i, j);
            }

            if self.should_report() {
                row_count += 1;
                if row_count % LOG_THRESHOLD == 0 { 
                    let c = self.pivots.count();
                    trace!("    [{row_count}/{total_rows}], {c} pivots.");
                }
            }
        }
     }

     #[cfg(feature = "multithread")]
     fn find_cycle_free_pivots_m(&mut self) {
        use yui::util::sync::SyncCounter;

        let remain_rows = self.remain_rows().collect_vec();
        let total_rows = remain_rows.len();

        trace!("  start find-cycle-free-pivots: {total_rows} rows");
        
        let n = self.cols();
        let pivots = RwLock::new(
            std::mem::take(&mut self.pivots)
        );
        let loc_pivots_tls = ThreadLocal::new();
        let loc_worker_tls = ThreadLocal::new();

        let report = self.should_report();
        let row_counter = SyncCounter::new();

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

            if report { 
                let row_count = row_counter.incr();            
                if row_count % LOG_THRESHOLD == 0 { 
                    let c = loc_pivots.count();
                    trace!("    [{row_count}/{total_rows}], {c} pivots.");
                }
            }
        });

        self.pivots = pivots.into_inner().unwrap();
     }

     #[cfg(feature = "multithread")]
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
        self.rows() > LOG_THRESHOLD && log::max_level() >= log::LevelFilter::Trace
     }
}

#[cfg(feature = "multithread")]
fn init_tls<T, F>(tl: &ThreadLocal<RefCell<T>>, f: F) -> &RefCell<T>
where T: Send, F: FnOnce() -> T {
    tl.get_or(|| RefCell::new( f() ) )
}

struct MatrixStr { 
    shape: (usize, usize),
    entries: Vec<Vec<Col>>,     // [row -> [col]]
    cands: Vec<AHashSet<Col>>,  // [row -> [col]]
    row_wght: Vec<f64>,         // [row -> weight]
    col_wght: Vec<f64>,         // [col -> weight]
}

impl MatrixStr { 
    fn new<R>(a: &SpMat<R>, piv_type: PivotType, pivot_cond: PivotCondition) -> Self
    where R: Ring, for<'x> &'x R: RingOps<R> { 
        let shape = match piv_type {
            PivotType::Rows => a.shape(),
            PivotType::Cols => (a.ncols(), a.nrows())
        };
        let t = match piv_type {
            PivotType::Rows => |i: usize, j: usize| (i, j),
            PivotType::Cols => |i, j| (j, i)
        };

        let (m, n) = shape;
        let mut entries = vec![vec![]; m];
        let mut row_wght = vec![0.0; m];
        let mut col_wght = vec![0.0; n];
        let mut cands = vec![AHashSet::new(); m];

        for (i, j, r) in a.iter() { 
            if r.is_zero() { continue }
            
            let (i, j) = t(i, j);
            entries[i].push(j);

            let w = r.c_weight();
            row_wght[i] += w;
            col_wght[j] += w;

            if pivot_cond.is_cand(r) { 
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
            a.ncols()
        } else { 
            a.nrows()
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

    #[allow(unused)]
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

#[repr(u8)]
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum EntryStatus { 
    None, Candidate, Occupied
}

struct RowWorker { 
    row: usize,
    status: Vec<EntryStatus>,
    ncand: usize,
    queue: VecDeque<Col>,
    queued: AHashSet<Col>
}

impl RowWorker {
    fn new(size: usize) -> Self { 
        let status = vec![EntryStatus::None; size];
        let queue = VecDeque::new();
        let queued = AHashSet::new();
        RowWorker {row: 0, status, ncand: 0, queue, queued }
    }

    fn clear(&mut self) {
        self.row = 0;
        self.status.fill(EntryStatus::None);
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

    #[allow(dead_code)]
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
        self.status[i] == EntryStatus::Candidate
    }

    fn set_candidate(&mut self, i: usize) { 
        assert_eq!(self.status[i], EntryStatus::None);
        self.status[i] = EntryStatus::Candidate;
        self.ncand += 1;
    }

    fn is_occupied(&self, i: usize) -> bool { 
        self.status[i] == EntryStatus::Occupied
    }

    fn set_occupied(&mut self, i: usize) { 
        if self.is_candidate(i) { 
            self.ncand -= 1;
        }
        self.status[i] = EntryStatus::Occupied;
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
        let a = SpMat::from_dense_data((6, 9), [
            1, 0, 1, 0, 0, 1, 1, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 2, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 1, 0, 3, 0, 0, 0, 0,
            0, 1, 0, 1, 0, 0, 1, 0, 1,
            1, 0, 1, 0, 1, 1, 0, 1, 1
        ]);
        let str = MatrixStr::new(&a, PivotType::Rows, PivotCondition::One);

        assert_eq!(str.entries, vec![
            vec![0,2,5,6,8], 
            vec![1,2,3,5,7], 
            vec![2,3,7,8], 
            vec![1,2,4], 
            vec![1,3,6,8], 
            vec![0,2,4,5,7,8]]
        );
        assert_eq!(str.row_wght, vec![5.0, 6.0, 4.0, 5.0, 4.0, 6.0]);
        assert_eq!(str.col_wght, vec![2.0, 3.0, 5.0, 3.0, 4.0, 3.0, 2.0, 4.0, 4.0]);
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
        let a = SpMat::from_dense_data((4, 4), [
            1, 0, 1, 0,
            0, 1, 1, 1,
            0, 0, 0, 0,
            0, 0, 1, 1,
        ]);
        let str = MatrixStr::new(&a, PivotType::Rows, PivotCondition::One);

        assert_eq!(str.head_col_in(0), Some(0));
        assert_eq!(str.head_col_in(1), Some(1));
        assert_eq!(str.head_col_in(2), None);
        assert_eq!(str.head_col_in(3), Some(2));
    }

    #[test]
    fn rows_cols() {
        let a = SpMat::<i32>::from_dense_data((4, 3), []);
        let pf = PivotFinder::new(&a, PivotType::Rows, PivotCondition::One);
        assert_eq!(pf.rows(), 4);
        assert_eq!(pf.cols(), 3);
    }

    #[test]
    fn pivot_data() { 
        let a = SpMat::from_dense_data((2, 4), [
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
        let a = SpMat::from_dense_data((4, 4), [
            1, 0, 1, 0,
            0, 1, 1, 1,
            0, 0, 0, 0,
            0, 0, 1, 1,
        ]);
        let mut pf = PivotFinder::new(&a, PivotType::Rows, PivotCondition::One);

        assert_eq!(pf.remain_rows().collect_vec(), vec![0,3,1]);

        pf.pivots.set(0, 0);

        assert_eq!(pf.remain_rows().collect_vec(), vec![3,1]);

        pf.pivots.set(1, 1);

        assert_eq!(pf.remain_rows().collect_vec(), vec![3]);
    }

    #[test]
    fn pivots() {
        let a = SpMat::from_dense_data((4, 4), [
            1, 0, 1, 0,
            0, 1, 1, 1,
            0, 0, 0, 0,
            0, 0, 1, 1,
        ]);
        let mut pf = PivotFinder::new(&a, PivotType::Rows, PivotCondition::One);

        assert_eq!(pf.occupied_cols(), AHashSet::new());

        pf.pivots.set(0, 0);

        assert_eq!(pf.occupied_cols(), AHashSet::from_iter([0,2]));

        pf.pivots.set(1, 1);

        assert_eq!(pf.occupied_cols(), AHashSet::from_iter([0,1,2,3]));
    }

    #[test]
    fn find_fl_pivots() {
        let a = SpMat::from_dense_data((6, 9), [
            1, 0, 1, 0, 0, 1, 1, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 1, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 1, 1
        ]);
        let mut pf = PivotFinder::new(&a, PivotType::Rows, PivotCondition::One);

        pf.find_fl_pivots();

        assert_eq!(pf.pivots.iter().collect_vec(), vec![(4, 2), (3, 1), (5, 5), (0, 0)]);
    }

    #[test]
    fn find_fl_col_pivots() { 
        let a = SpMat::from_dense_data((6, 9), [
            1, 0, 0, 0, 0, 1, 0, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 1, 0, 1, 0
        ]);
        let mut pf = PivotFinder::new(&a, PivotType::Rows, PivotCondition::One);

        pf.find_fl_col_pivots();

        assert_eq!(pf.pivots.iter().collect_vec(), vec![(4, 2), (3, 4), (0, 0), (5, 7), (2, 3)]);
    }

    #[test]
    fn find_fl_row_col_pivots() { 
        let a = SpMat::from_dense_data((6, 9), [
            1, 0, 0, 0, 0, 1, 0, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 1, 0, 1, 0
        ]);
        let mut pf = PivotFinder::new(&a, PivotType::Rows, PivotCondition::One);

        pf.find_fl_pivots();

        assert_eq!(pf.pivots.iter().collect_vec(), vec![(4, 2), (3, 1), (0, 0)]);

        pf.find_fl_col_pivots();

        assert_eq!(pf.pivots.iter().collect_vec(), vec![(4, 2), (3, 1), (0, 0), (5, 7), (2, 3)]);
    }

    #[test]
    fn find_cycle_free_pivots_s() {
        let a = SpMat::from_dense_data((6, 9), [
            1, 0, 0, 0, 0, 1, 0, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 1, 0, 1, 0
        ]);
        let mut pf = PivotFinder::new(&a, PivotType::Rows, PivotCondition::One);

        pf.find_cycle_free_pivots_s();

        assert_eq!(pf.pivots.iter().collect_vec(), vec![(4, 2), (3, 4), (0, 0), (5, 1), (2, 3)]);
    }

    #[cfg(feature = "multithread")]
    #[test]
    fn find_cycle_free_pivots_m() {
        let a = SpMat::from_dense_data((6, 9), [
            1, 0, 0, 0, 0, 1, 0, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 1, 0, 1, 0
        ]);
        let mut pf = PivotFinder::new(&a, PivotType::Rows, PivotCondition::One);

        pf.find_cycle_free_pivots_m();

        assert!(pf.pivots.count() >= 5);
    }

    #[test]
    fn zero() { 
        let a = SpMat::from_dense_data((1, 1), [0]);
        let pivs = find_pivots(&a, PivotType::Rows, PivotCondition::One);
        let r = pivs.len();
        assert_eq!(r, 0);
    }

    #[test]
    fn id_1() { 
        let a = SpMat::from_dense_data((1, 1), [1]);
        let pivs = find_pivots(&a, PivotType::Rows, PivotCondition::One);
        let r = pivs.len();
        assert_eq!(r, 1);
    }

    #[test]
    fn id_2() { 
        let a = SpMat::from_dense_data((2, 2), [
            1, 0, 0, 1
        ]);
        let pivs = find_pivots(&a, PivotType::Rows, PivotCondition::One);
        let r = pivs.len();
        assert_eq!(r, 2);
    }

    #[test]
    fn result() { 
        let a = SpMat::from_dense_data((6, 9), [
            1, 0, 0, 0, 0, 1, 0, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 1, 0, 1, 0
        ]);
        let pivs = find_pivots(&a, PivotType::Rows, PivotCondition::One);
        let r = pivs.len();
        assert_eq!(r, 5);
        
        let (p, q) = perms_by_pivots(&a, &pivs);
        let b = a.permute(p.view(), q.view()).into_dense();

        assert!((0..r).all(|i| b[(i, i)].is_one()));
        assert!((0..r).all(|j| {
            (j+1..r).all(|i| b[(i, j)].is_zero())
        }));
    }

    #[test]
    fn result_cols() { 
        let a = SpMat::from_dense_data((6, 9), [
            1, 0, 0, 0, 0, 1, 0, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 1, 0, 1, 0
        ]);
        let pivs = find_pivots(&a, PivotType::Cols, PivotCondition::One);
        let r = pivs.len();
        assert_eq!(r, 6);
        
        let (p, q) = perms_by_pivots(&a, &pivs);
        let b = a.permute(p.view(), q.view()).into_dense();

        assert!((0..r).all(|i| b[(i, i)].is_one()));
        assert!((0..r).all(|i| {
            (i+1..r).all(|j| b[(i, j)].is_zero())
        }));
    }

    #[test]
    fn rand() {
        let d = 0.1;
        let shape = (60, 80);
        let a = SpMat::<i32>::rand(shape, d);

        let pivs = find_pivots(&a, PivotType::Rows, PivotCondition::One);
        let r = pivs.len();
        assert!(r > 10);
        
        let (p, q) = perms_by_pivots(&a, &pivs);
        let b = a.permute(p.view(), q.view()).into_dense();

        assert!((0..r).all(|i| b[(i, i)].is_one()));
        assert!((0..r).all(|j| {
            (j+1..r).all(|i| b[(i, j)].is_zero())
        }))
    }
}