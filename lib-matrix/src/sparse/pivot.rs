//! Heuristic pivot finder for sparse PLUQ.
//!
//! Implementation based on:
//!
//! - "Parallel Sparse PLUQ Factorization modulo p", Charles Bouillaguet,
//!   Claire Delaplace, Marie-Emilie Voge.
//!   <https://hal.inria.fr/hal-01646133/document>
//! - see also: SpaSM (Sparse direct Solver Modulo p),
//!   <https://github.com/cbouilla/spasm>.

use std::cmp::Ordering;
use std::collections::VecDeque;
use itertools::Itertools;
use log::*;

use yui_core::{Ring, RingOps};
use yui_core::algo::TopSort;
use crate::Perm;
use super::*;

cfg_if::cfg_if! {
    if #[cfg(feature = "multithread")] {
        use std::cell::RefCell;
        use std::sync::RwLock;
        use thread_local::ThreadLocal;
        use rayon::prelude::*;
        use yui_core::util::sync::SyncCounter;
    }
}

const LOG_THRESHOLD: usize = 10_000;
const DEFAULT_MAX_PIVOT: usize = usize::MAX;

/// Whether pivots are picked along rows (each pivot eliminates a row's
/// other entries) or along columns. Affects how the resulting `L`/`U`
/// blocks are oriented after reduction.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum PivotType {
    Rows, Cols
}

impl PivotType {
    fn str(&self) -> &'static str { 
        match self {
            PivotType::Rows => "row",
            PivotType::Cols => "col"
        }
    }
}

/// Which entries are eligible as pivots:
/// - `One` — only `±1`.
/// - `Weight(w)` — any unit whose `c_weight()` is at most `w`.
/// - `AnyUnit` — any unit of the ring.
#[derive(Clone, Copy, Debug)]
pub enum PivotCondition {
    One, Weight(f64), AnyUnit
}

/// Knobs for [`find_pivots`].
#[derive(Clone, Copy, Debug)]
pub struct PivotFinderConfig {
    pub piv_type: PivotType,
    pub piv_cond: PivotCondition,
    pub max_pivots: usize,
}

impl Default for PivotFinderConfig {
    fn default() -> Self {
        Self {
            piv_type: PivotType::Rows,
            piv_cond: PivotCondition::One,
            max_pivots: DEFAULT_MAX_PIVOT,
        }
    }
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

/// Searches for pivots in `a` and returns row/column permutations `(p, q)`
/// that bring the chosen pivots to the leading r×r block of `p·a·q⁻¹`, along
/// with the rank `r` (number of pivots found).
pub fn find_pivots<R>(a: &SpMat<R>, config: PivotFinderConfig) -> (Perm, Perm, usize)
where R: Ring, for<'x> &'x R: RingOps<R> {
    let (m, n) = a.shape();

    if a.is_zero() {
        return (Perm::id(m), Perm::id(n), 0);
    }

    debug!("find {} pivots: {:?}", config.piv_type.str(), a.shape());

    let mut pf = PivotFinder::new(a, &config);
    pf.find_pivots();
    let pivs = pf.result();

    debug!("  found {} {} pivots", pivs.len(), config.piv_type.str());

    let p = Perm::forward_indices(m, pivs.iter().map(|(i, _)| *i));
    let q = Perm::forward_indices(n, pivs.iter().map(|(_, j)| *j));
    let r = pivs.len();
    
    (p, q, r)
}

type Row = usize;
type Col = usize;

pub struct PivotFinder {
    str: MatrixStr,
    pivots: PivotData,
    piv_type: PivotType,
    max_pivots: usize,
}

impl PivotFinder {
    pub fn new<R>(a: &SpMat<R>, config: &PivotFinderConfig) -> Self
    where R: Ring, for<'x> &'x R: RingOps<R> {
        let str = MatrixStr::new(a, config.piv_type, config.piv_cond);
        let pivots = PivotData::new(a, config.piv_type);
        PivotFinder { str, pivots, piv_type: config.piv_type, max_pivots: config.max_pivots }
    }

    pub fn find_pivots(&mut self) {
        trace!("pivots: {:?} ..", self.str.shape());

        self.find_fl_pivots();
        if self.pivots.count() < self.max_pivots {
            self.find_fl_col_pivots();
        }
        if self.pivots.count() < self.max_pivots {
            self.find_cycle_free_pivots();
        }

        trace!("pivots: {:?} => {}.", self.str.shape(), self.pivots.count());
    }

    pub fn result(&self) -> Vec<(usize, usize)> {
        let mut ts = TopSort::new();
        for (i, j) in self.pivots.iter() {
            ts.add_node(j);
            for j2 in self.str.cols_in(i).filter(|&j2| j != j2 && self.pivots.has_col(j2)) {
                ts.add_edge(j, j2);
            }
        }
        let sorted = ts.into_sorted().unwrap();
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

    fn remain_rows(&self) -> impl Iterator<Item = Row> + '_ {
        self.str.target_rows.iter().copied()
            .filter(|&i| !self.pivots.is_piv_row(i))
    }

    // occupancy bitmap over columns: the finder already holds O(n_cols) state (col_wght, RowWorker),
    // so the bitmap is free, and lookups are one per matrix entry — they must be cheaper than hashing.
    fn occupied_cols(&self) -> Vec<bool> {
        let mut occ = vec![false; self.cols()];
        for (i, _) in self.pivots.iter() {
            for j in self.str.cols_in(i) {
                occ[j] = true;
            }
        }
        occ
    }

    fn find_fl_pivots(&mut self) {
        let remain_rows: Vec<_> = self.remain_rows().collect();

        for i in remain_rows {
            let Some((j, is_cand)) = self.str.head(i) else { continue };

            if is_cand && !self.pivots.has_col(j) {
                self.pivots.set(i, j);
                if self.pivots.count() >= self.max_pivots { break; }
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

            for (j, is_cand) in self.str.entries_in(i) {
                if is_cand && !occ_cols[j] {
                    cands.push(j);
                }
            }

            let Some(j) = cands.into_iter().min_by(|&j1, &j2|
                self.str.cmp_cols(j1, j2)
            ) else { continue };

            self.pivots.set(i, j);

            for j in self.str.cols_in(i) {
                occ_cols[j] = true;
            }

            if self.pivots.count() >= self.max_pivots { break; }
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
            if self.pivots.count() >= self.max_pivots { break; }

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
        let remain_rows = self.remain_rows().collect_vec();
        let total_rows = remain_rows.len();

        trace!("  start find-cycle-free-pivots: {total_rows} rows");

        let n = self.cols();
        let count = SyncCounter::new(self.pivots.count()); // lock-free mirror of `pivots.count()`
        let pivots = RwLock::new(
            std::mem::take(&mut self.pivots)
        );
        let loc_pivots_tls = ThreadLocal::new();
        let loc_worker_tls = ThreadLocal::new();

        let report = self.should_report();
        let row_counter = SyncCounter::new(0);

        remain_rows.par_iter().for_each(|&i| {
            if count.count() >= self.max_pivots { return; }

            let mut loc_pivots = init_tls(&loc_pivots_tls, ||
                pivots.read().unwrap().clone()
            ).borrow_mut();

            let mut w = init_tls(&loc_worker_tls, ||
                RowWorker::new(n)
            ).borrow_mut();

            loc_pivots.update_from(&pivots.read().unwrap());
            w.init(i, &self.str, &loc_pivots);

            self.find_cycle_free_pivots_in(&pivots, &count, &mut loc_pivots, &mut w);

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
     fn find_cycle_free_pivots_in(&self, pivots: &RwLock<PivotData>, count: &SyncCounter, loc_pivots: &mut PivotData, w: &mut RowWorker) {
        loop {
            w.traverse(&self.str, loc_pivots);

            let Some(j) = w.choose_candidate(&self.str) else {
                break
            };

            // If changes are made in other threads, update `loc_pivots` and retry.
            // Otherwise, modify `pivots` and exit.

            let mut pivots = pivots.write().unwrap();
            w.update_diff(loc_pivots, &pivots);

            if w.should_retry() {
                loc_pivots.update_from(&pivots);
                continue
            } else {
                if pivots.count() < self.max_pivots {
                    pivots.set(w.row, j);
                    count.set(pivots.count());
                }
                break
            }
        }
     }

     fn should_report(&self) -> bool { 
        self.rows() > LOG_THRESHOLD && log::max_level() >= log::LevelFilter::Debug
     }
}

#[cfg(feature = "multithread")]
fn init_tls<T, F>(tl: &ThreadLocal<RefCell<T>>, f: F) -> &RefCell<T>
where T: Send, F: FnOnce() -> T {
    tl.get_or(|| RefCell::new( f() ) )
}

struct MatrixStr {
    shape: (usize, usize),
    entries: Vec<Vec<(Col, bool)>>,    // [row -> sorted [(col, is_cand)]]
    target_rows: Vec<Row>,             // non-empty rows, sorted by (row_wght, row index)
    col_wght: Vec<f64>,                // [col -> weight]
}

impl MatrixStr {
    fn new<R>(a: &SpMat<R>, piv_type: PivotType, pivot_cond: PivotCondition) -> Self
    where R: Ring, for<'x> &'x R: RingOps<R> {
        let shape = match piv_type {
            PivotType::Rows => a.shape(),
            PivotType::Cols => (a.n_cols(), a.n_rows())
        };
        let t = match piv_type {
            PivotType::Rows => |i: usize, j: usize| (i, j),
            PivotType::Cols => |i, j| (j, i)
        };

        let (m, n) = shape;
        let mut entries = vec![vec![]; m];
        let mut row_wght = vec![0.0; m];
        let mut col_wght = vec![0.0; n];

        for (i, j, r) in a.iter() {
            if r.is_zero() { continue }

            let (i, j) = t(i, j);
            entries[i].push((j, pivot_cond.is_cand(r)));

            let w = r.c_weight();
            row_wght[i] += w;
            col_wght[j] += w;
        }

        let mut target_rows: Vec<Row> = (0..m).filter(|&i| !entries[i].is_empty()).collect();
        target_rows.sort_unstable_by(|&i1, &i2| {
            row_wght[i1].partial_cmp(&row_wght[i2])
                .unwrap_or(Ordering::Equal)
                .then(i1.cmp(&i2))
        });

        Self { shape, entries, col_wght, target_rows }
    }

    fn shape(&self) -> (usize, usize) {
        self.shape
    }

    fn head(&self, i: Row) -> Option<(Col, bool)> {
        self.entries[i].first().copied()
    }

    fn entries_in(&self, i: Row) -> impl Iterator<Item = (Col, bool)> + '_ {
        self.entries[i].iter().copied()
    }

    fn cols_in(&self, i: Row) -> impl Iterator<Item = Col> + '_ {
        self.entries[i].iter().map(|(c, _)| *c)
    }

    fn cmp_cols(&self, j1: Col, j2: Col) -> Ordering {
        if let Some(o) = self.col_wght[j1].partial_cmp(&self.col_wght[j2]) {
            o.then(Ord::cmp(&j1, &j2))
        } else {
            Ordering::Equal
        }
    }
}

#[derive(Clone, Default)]
struct PivotData {
    data: Vec<Option<Row>>,   // col -> row
    indices: Vec<Col>,
    is_piv_row: Vec<bool>,    // row -> is this a pivot row?
}

impl PivotData {
    fn new<R>(a: &SpMat<R>, piv_type: PivotType) -> Self
    where R: Ring, for<'x> &'x R: RingOps<R> {
        let (m, n) = match piv_type {
            PivotType::Rows => (a.n_rows(), a.n_cols()),
            PivotType::Cols => (a.n_cols(), a.n_rows()),
        };
        let data = vec![None; n];
        let indices = vec![];
        let is_piv_row = vec![false; m];
        Self { data, indices, is_piv_row }
    }

    fn count(&self) -> usize {
        self.indices.len()
    }

    fn has_col(&self, j: Col) -> bool {
        self.data[j].is_some()
    }

    fn is_piv_row(&self, i: Row) -> bool {
        self.is_piv_row[i]
    }

    fn row_for(&self, j: Col) -> Option<Row> {
        self.data[j]
    }

    fn set(&mut self, i: Row, j: Col) {
        assert!(!self.has_col(j));
        self.data[j] = Some(i);
        self.indices.push(j);
        self.is_piv_row[i] = true;
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
    candidates: Vec<Col>,    // columns ever marked Candidate (may include stale entries reclassified to Occupied)
    queue: VecDeque<Col>,
    queued: Vec<bool>,       // queued[j] = true iff j has ever been pushed to `queue`
    touched: Vec<Col>,       // columns whose status/queued changed — `clear` resets only these, O(touched) not O(n)
}

impl RowWorker {
    fn new(size: usize) -> Self {
        let status = vec![EntryStatus::None; size];
        let candidates = Vec::new();
        let queue = VecDeque::new();
        let queued = vec![false; size];
        let touched = Vec::new();
        RowWorker { row: 0, status, ncand: 0, candidates, queue, queued, touched }
    }

    fn clear(&mut self) {
        self.row = 0;
        self.ncand = 0;
        for &j in self.touched.iter() {
            self.status[j] = EntryStatus::None;
            self.queued[j] = false;
        }
        self.touched.clear();
        self.candidates.clear();
        self.queue.clear();
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

        for (j, is_cand) in str.entries_in(i) {
            if pivots.has_col(j) {
                self.enqueue(j);
                self.set_occupied(j);
            } else if is_cand {
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

            for j2 in str.cols_in(i2) {
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
        if self.ncand == 0 { return None; }
        self.candidates.iter().copied()
            .filter(|&j| self.is_candidate(j))
            .min_by(|&j1, &j2|
                str.cmp_cols(j1, j2)
            )
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
        self.candidates.push(i);
        self.touched.push(i);
        self.ncand += 1;
    }

    fn is_occupied(&self, i: usize) -> bool {
        self.status[i] == EntryStatus::Occupied
    }

    fn set_occupied(&mut self, i: usize) {
        match self.status[i] {
            EntryStatus::Candidate => {
                self.ncand -= 1; // already in `touched` via set_candidate
            }
            EntryStatus::None => {
                self.touched.push(i);
            }
            EntryStatus::Occupied => {
                return;
            }
        }
        self.status[i] = EntryStatus::Occupied;
    }

    fn enqueue(&mut self, i: Col) {
        self.queue.push_back(i);
        if !self.queued[i] {
            self.queued[i] = true;
            self.touched.push(i);
        }
    }

    fn dequeue(&mut self) -> Option<Col> {
        self.queue.pop_front()
    }

    fn is_queued(&self, i: Col) -> bool {
        self.queued[i]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::{Zero, One};
 
    #[test]
    fn str_init() {
        let a = SpMat::from_row_major((6, 9), [
            1, 0, 1, 0, 0, 1, 1, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 2, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 1, 0, 3, 0, 0, 0, 0,
            0, 1, 0, 1, 0, 0, 1, 0, 1,
            1, 0, 1, 0, 1, 1, 0, 1, 1
        ]);
        let str = MatrixStr::new(&a, PivotType::Rows, PivotCondition::One);

        assert_eq!(str.entries, vec![
            vec![(0,true), (2,true), (5,true), (6,true), (8,true)],
            vec![(1,true), (2,true), (3,true), (5,true), (7,false)],
            vec![(2,true), (3,true), (7,true), (8,true)],
            vec![(1,true), (2,true), (4,false)],
            vec![(1,true), (3,true), (6,true), (8,true)],
            vec![(0,true), (2,true), (4,true), (5,true), (7,true), (8,true)],
        ]);
        assert_eq!(str.col_wght, vec![2.0, 3.0, 5.0, 3.0, 4.0, 3.0, 2.0, 4.0, 4.0]);
        assert_eq!(str.target_rows, vec![2, 4, 0, 3, 1, 5]);
    }

    #[test]
    fn str_row_head() {
        let a = SpMat::from_row_major((4, 4), [
            1, 0, 1, 0,
            0, 1, 1, 1,
            0, 0, 0, 0,
            0, 0, 1, 1,
        ]);
        let str = MatrixStr::new(&a, PivotType::Rows, PivotCondition::One);

        assert_eq!(str.head(0), Some((0, true)));
        assert_eq!(str.head(1), Some((1, true)));
        assert_eq!(str.head(2), None);
        assert_eq!(str.head(3), Some((2, true)));
    }

    #[test]
    fn rows_cols() {
        let a = SpMat::<i32>::from_row_major((4, 3), []);
        let pf = PivotFinder::new(&a, &Default::default());
        assert_eq!(pf.rows(), 4);
        assert_eq!(pf.cols(), 3);
    }

    #[test]
    fn pivot_data() { 
        let a = SpMat::from_row_major((2, 4), [
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
        let a = SpMat::from_row_major((4, 4), [
            1, 0, 1, 0,
            0, 1, 1, 1,
            0, 0, 0, 0,
            0, 0, 1, 1,
        ]);
        let mut pf = PivotFinder::new(&a, &Default::default());

        assert_eq!(pf.remain_rows().collect_vec(), vec![0,3,1]);

        pf.pivots.set(0, 0);

        assert_eq!(pf.remain_rows().collect_vec(), vec![3,1]);

        pf.pivots.set(1, 1);

        assert_eq!(pf.remain_rows().collect_vec(), vec![3]);
    }

    #[test]
    fn pivots() {
        let a = SpMat::from_row_major((4, 4), [
            1, 0, 1, 0,
            0, 1, 1, 1,
            0, 0, 0, 0,
            0, 0, 1, 1,
        ]);
        let mut pf = PivotFinder::new(&a, &Default::default());

        let occ = |pf: &PivotFinder| pf.occupied_cols().iter().positions(|&b| b).collect_vec();

        assert_eq!(occ(&pf), vec![]);

        pf.pivots.set(0, 0);

        assert_eq!(occ(&pf), vec![0, 2]);

        pf.pivots.set(1, 1);

        assert_eq!(occ(&pf), vec![0, 1, 2, 3]);
    }

    #[test]
    fn find_fl_pivots() {
        let a = SpMat::from_row_major((6, 9), [
            1, 0, 1, 0, 0, 1, 1, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 1, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 1, 1
        ]);
        let mut pf = PivotFinder::new(&a, &Default::default());

        pf.find_fl_pivots();

        assert_eq!(pf.pivots.iter().collect_vec(), vec![(4, 2), (3, 1), (5, 5), (0, 0)]);
    }

    #[test]
    fn find_fl_col_pivots() { 
        let a = SpMat::from_row_major((6, 9), [
            1, 0, 0, 0, 0, 1, 0, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 1, 0, 1, 0
        ]);
        let mut pf = PivotFinder::new(&a, &Default::default());

        pf.find_fl_col_pivots();

        assert_eq!(pf.pivots.iter().collect_vec(), vec![(4, 2), (3, 4), (0, 0), (5, 7), (2, 3)]);
    }

    #[test]
    fn find_fl_row_col_pivots() { 
        let a = SpMat::from_row_major((6, 9), [
            1, 0, 0, 0, 0, 1, 0, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 1, 0, 1, 0
        ]);
        let mut pf = PivotFinder::new(&a, &Default::default());

        pf.find_fl_pivots();

        assert_eq!(pf.pivots.iter().collect_vec(), vec![(4, 2), (3, 1), (0, 0)]);

        pf.find_fl_col_pivots();

        assert_eq!(pf.pivots.iter().collect_vec(), vec![(4, 2), (3, 1), (0, 0), (5, 7), (2, 3)]);
    }

    #[test]
    fn find_cycle_free_pivots_s() {
        let a = SpMat::from_row_major((6, 9), [
            1, 0, 0, 0, 0, 1, 0, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 1, 0, 1, 0
        ]);
        let mut pf = PivotFinder::new(&a, &Default::default());

        pf.find_cycle_free_pivots_s();

        assert_eq!(pf.pivots.iter().collect_vec(), vec![(4, 2), (3, 4), (0, 0), (5, 1), (2, 3)]);
    }

    #[cfg(feature = "multithread")]
    #[test]
    fn find_cycle_free_pivots_m() {
        let a = SpMat::from_row_major((6, 9), [
            1, 0, 0, 0, 0, 1, 0, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 1, 0, 1, 0
        ]);
        let mut pf = PivotFinder::new(&a, &Default::default());

        pf.find_cycle_free_pivots_m();

        assert!(pf.pivots.count() >= 5);
    }

    #[test]
    fn zero() {
        let a = SpMat::from_row_major((1, 1), [0]);
        let (p, q, r) = find_pivots(&a, Default::default());
        assert_eq!(r, 0);
        assert_eq!(p.len(), 1);
        assert_eq!(q.len(), 1);
        assert!(p.is_id());
        assert!(q.is_id());
    }

    #[test]
    fn id_1() {
        let a = SpMat::from_row_major((1, 1), [1]);
        let (p, q, r) = find_pivots(&a, Default::default());
        assert_eq!(r, 1);
        assert_eq!(p.len(), 1);
        assert_eq!(q.len(), 1);
        assert!(p.is_id());
        assert!(q.is_id());
    }

    #[test]
    fn id_2() {
        let a = SpMat::from_row_major((2, 2), [
            1, 0, 0, 1
        ]);
        let (p, q, r) = find_pivots(&a, Default::default());
        assert_eq!(r, 2);
        assert_eq!(p.len(), 2);
        assert_eq!(q.len(), 2);
    }

    #[test]
    fn result() {
        let a = SpMat::from_row_major((6, 9), [
            1, 0, 0, 0, 0, 1, 0, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 1, 0, 1, 0
        ]);
        let (p, q, r) = find_pivots(&a, Default::default());
        assert_eq!(r, 5);
        assert_eq!(p.len(), 6);
        assert_eq!(q.len(), 9);

        let b = a.permute(&p, &q).into_dense();

        assert!((0..r).all(|i| b[(i, i)].is_one()));
        assert!((0..r).all(|j| {
            (j+1..r).all(|i| b[(i, j)].is_zero())
        }));
    }

    #[test]
    fn result_cols() {
        let a = SpMat::from_row_major((6, 9), [
            1, 0, 0, 0, 0, 1, 0, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 1, 0, 1, 0
        ]);
        let config = PivotFinderConfig { piv_type: PivotType::Cols, ..Default::default() };
        let (p, q, r) = find_pivots(&a, config);
        assert_eq!(r, 6);
        assert_eq!(p.len(), 6);
        assert_eq!(q.len(), 9);

        let b = a.permute(&p, &q).into_dense();

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

        let (p, q, r) = find_pivots(&a, Default::default());
        assert!(r > 10);
        assert_eq!(p.len(), shape.0);
        assert_eq!(q.len(), shape.1);

        let b = a.permute(&p, &q).into_dense();

        assert!((0..r).all(|i| b[(i, i)].is_one()));
        assert!((0..r).all(|j| {
            (j+1..r).all(|i| b[(i, j)].is_zero())
        }))
    }

    #[test]
    fn max_pivots() {
        // Full rank of this matrix is 5; limiting to 3 must return ≤ 3 pivots.
        let a = SpMat::from_row_major((6, 9), [
            1, 0, 0, 0, 0, 1, 0, 0, 1,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 0, 0, 1, 1,
            0, 1, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 1, 0, 1, 0
        ]);
        let config = PivotFinderConfig { max_pivots: 3, ..Default::default() };
        let (p, q, r) = find_pivots(&a, config);
        assert!(r <= 3);
        assert_eq!(p.len(), 6);
        assert_eq!(q.len(), 9);

        let b = a.permute(&p, &q).into_dense();
        assert!((0..r).all(|i| b[(i, i)].is_one()));
        assert!((0..r).all(|j| {
            (j+1..r).all(|i| b[(i, j)].is_zero())
        }));
    }

}