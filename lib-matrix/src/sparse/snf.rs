//! Sparse Smith normal form over a Euclidean ring.
//!
//! Mirrors [`crate::dense::snf`] (same pivoting and Euclidean elimination),
//! but operates on hashed sparse rows with a column index, so large
//! reducer-residuals never materialize as dense matrices. The four
//! transformation matrices are reconstructed at the end by replaying the
//! recorded elementary operations on sparse identities.

use log::debug;
use rustc_hash::{FxHashMap, FxHashSet};
use yui_core::{EucRing, EucRingOps};

use crate::MatTrait;
use crate::dense::snf::SnfFlags;
use super::SpMat;

/// Result of a sparse SNF: `p * a * q = result` (diagonal), `pinv * result * qinv = a`.
pub struct SpSnf<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    result: SpMat<R>,
    diag: Vec<R>,
    p:    Option<SpMat<R>>,
    pinv: Option<SpMat<R>>,
    q:    Option<SpMat<R>>,
    qinv: Option<SpMat<R>>,
}

impl<R> SpSnf<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn result(&self) -> &SpMat<R> {
        &self.result
    }

    pub fn p(&self) -> Option<&SpMat<R>> {
        self.p.as_ref()
    }

    pub fn pinv(&self) -> Option<&SpMat<R>> {
        self.pinv.as_ref()
    }

    pub fn q(&self) -> Option<&SpMat<R>> {
        self.q.as_ref()
    }

    pub fn qinv(&self) -> Option<&SpMat<R>> {
        self.qinv.as_ref()
    }

    pub fn rank(&self) -> usize {
        self.diag.len()
    }

    pub fn factors(&self) -> Vec<&R> {
        self.diag.iter().collect()
    }
}

/// Below this core area (non-zero rows × non-zero cols), the matrix is
/// compacted by permutation and passed to the dense SNF instead.
const DENSE_SNF_MAX_AREA: usize = 65536; // 256 × 256

/// Computes the sparse SNF of `a`, producing the transformation matrices
/// selected by `flags = [p, pinv, q, qinv]`.
///
/// Small inputs take a fast pass: permute the non-zero core to the front and
/// run [`crate::dense::snf`] on it; large ones are eliminated sparsely.
pub fn sp_snf<R>(a: &SpMat<R>, flags: SnfFlags) -> SpSnf<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    sp_snf_with(a, flags, DENSE_SNF_MAX_AREA)
}

fn sp_snf_with<R>(a: &SpMat<R>, flags: SnfFlags, dense_max_area: usize) -> SpSnf<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    let (row_idx, col_idx) = nz_indices(a);
    let (m0, n0) = (row_idx.len(), col_idx.len());

    debug!("start sparse snf: {:?}, nnz: {}, core: {m0}x{n0}, flags: {:?}", a.shape(), a.iter_nz().count(), flags);

    if m0 * n0 <= dense_max_area {
        return dense_snf_in(a, flags, row_idx, col_idx);
    }

    let mut calc = SpSnfCalc::new(a, flags);
    calc.process();

    debug!("sparse snf done, rank: {}", calc.rank);

    calc.into_result()
}

// Sorted non-zero row / col indices of `a`.
fn nz_indices<R>(a: &SpMat<R>) -> (Vec<usize>, Vec<usize>)
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    let mut rows = vec![false; a.n_rows()];
    let mut cols = vec![false; a.n_cols()];
    for (i, j, _) in a.iter_nz() {
        rows[i] = true;
        cols[j] = true;
    }
    let row_idx = rows.iter().enumerate().filter_map(|(i, &b)| b.then_some(i)).collect();
    let col_idx = cols.iter().enumerate().filter_map(|(j, &b)| b.then_some(j)).collect();
    (row_idx, col_idx)
}

// Fast pass: gather the non-zero core `a[row_idx, col_idx]` to the front,
// run the dense SNF on it, and lift. With `R`/`C` the gathering permutations,
// `R·a·C = [[core, 0], [0, 0]]`, so
// `p = diag(p_d, I)·R`, `q = C·diag(q_d, I)`, and inverses transposed-fashion.
fn dense_snf_in<R>(a: &SpMat<R>, flags: SnfFlags, row_idx: Vec<usize>, col_idx: Vec<usize>) -> SpSnf<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    use crate::dense::Mat;
    use crate::dense::snf::snf_in_place;

    let (m, n) = a.shape();
    let (m0, n0) = (row_idx.len(), col_idx.len());

    debug!("  snf fast-pass: dense on {m0}x{n0} core");

    let row_pos: FxHashMap<usize, usize> = row_idx.iter().enumerate().map(|(k, &i)| (i, k)).collect();
    let col_pos: FxHashMap<usize, usize> = col_idx.iter().enumerate().map(|(k, &j)| (j, k)).collect();

    let mut core = Mat::zero((m0, n0));
    for (i, j, v) in a.iter_nz() {
        core[(row_pos[&i], col_pos[&j])] = v.clone();
    }

    let s = snf_in_place(core, flags);
    let rank = s.rank();
    let (res_d, [p_d, pinv_d, q_d, qinv_d]) = s.destruct();

    let diag: Vec<R> = (0..rank).map(|i| res_d[(i, i)].clone()).collect();
    let result = SpMat::from_entries((m, n), diag.iter().enumerate().map(|(i, v)| (i, i, v.clone())));

    let rest_rows = || (0..m).filter(|i| !row_pos.contains_key(i));
    let rest_cols = || (0..n).filter(|j| !col_pos.contains_key(j));

    // p[k, row_idx[l]] = p_d[k, l];  p[m0+t, rest_row_t] = 1
    let p = p_d.map(|p_d| SpMat::from_entries((m, m),
        dense_nz(&p_d).map(|(k, l, v)| (k, row_idx[l], v))
            .chain(rest_rows().enumerate().map(|(t, i)| (m0 + t, i, R::one())))
    ));
    // pinv[row_idx[k], l] = pinv_d[k, l];  pinv[rest_row_t, m0+t] = 1
    let pinv = pinv_d.map(|pinv_d| SpMat::from_entries((m, m),
        dense_nz(&pinv_d).map(|(k, l, v)| (row_idx[k], l, v))
            .chain(rest_rows().enumerate().map(|(t, i)| (i, m0 + t, R::one())))
    ));
    // q[col_idx[k], l] = q_d[k, l];  q[rest_col_t, n0+t] = 1
    let q = q_d.map(|q_d| SpMat::from_entries((n, n),
        dense_nz(&q_d).map(|(k, l, v)| (col_idx[k], l, v))
            .chain(rest_cols().enumerate().map(|(t, j)| (j, n0 + t, R::one())))
    ));
    // qinv[k, col_idx[l]] = qinv_d[k, l];  qinv[n0+t, rest_col_t] = 1
    let qinv = qinv_d.map(|qinv_d| SpMat::from_entries((n, n),
        dense_nz(&qinv_d).map(|(k, l, v)| (k, col_idx[l], v))
            .chain(rest_cols().enumerate().map(|(t, j)| (n0 + t, j, R::one())))
    ));

    SpSnf { result, diag, p, pinv, q, qinv }
}

fn dense_nz<R>(a: &crate::dense::Mat<R>) -> impl Iterator<Item = (usize, usize, R)> + '_
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    let (m, n) = a.shape();
    (0..m).flat_map(move |i| (0..n).filter_map(move |j| {
        let v = &a[(i, j)];
        (!v.is_zero()).then(|| (i, j, v.clone()))
    }))
}

// Recorded elementary operations, replayed to build the trans matrices.
// `Left([a,b,c,d], i, j)`: rows (i, j) ← (a·rᵢ + b·rⱼ, c·rᵢ + d·rⱼ), det = 1.
// `Right([a,b,c,d], i, j)`: cols (i, j) ← (a·cᵢ + b·cⱼ, c·cᵢ + d·cⱼ), det = 1.
enum Op<R> {
    SwapRows(usize, usize),
    SwapCols(usize, usize),
    MulRow(usize, R),
    MulCol(usize, R),
    Left([R; 4], usize, usize),
    Right([R; 4], usize, usize),
}

struct SpSnfCalc<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    shape: (usize, usize),
    rows: Vec<FxHashMap<usize, R>>,  // row -> (col -> value)
    cols: Vec<FxHashSet<usize>>,     // col -> rows with a non-zero entry
    ops: Vec<Op<R>>,
    flags: SnfFlags,
    rank: usize,
}

impl<R> SpSnfCalc<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    fn new(a: &SpMat<R>, flags: SnfFlags) -> Self {
        let (m, n) = a.shape();
        let mut rows = vec![FxHashMap::default(); m];
        let mut cols = vec![FxHashSet::default(); n];

        for (i, j, v) in a.iter_nz() {
            rows[i].insert(j, v.clone());
            cols[j].insert(i);
        }

        Self { shape: (m, n), rows, cols, ops: Vec::new(), flags, rank: 0 }
    }

    fn process(&mut self) {
        self.eliminate_all();
        debug!("  snf eliminate done: rank {}, ops {}; normalize..", self.rank, self.ops.len());
        self.diag_normalize();
    }

    fn into_result(self) -> SpSnf<R> {
        let diag: Vec<R> = (0..self.rank).map(|i| self.rows[i][&i].clone()).collect();
        let result = SpMat::from_entries(self.shape, self.rows.iter().enumerate().flat_map(|(i, row)|
            row.iter().map(move |(&j, v)| (i, j, v.clone()))
        ));
        let [p, pinv, q, qinv] = self.make_trans();
        SpSnf { result, diag, p, pinv, q, qinv }
    }

    // ---- elimination ----

    fn eliminate_all(&mut self) {
        let (m, n) = self.shape;
        let mut i = 0;

        for j in 0..n {
            if i >= m { break }
            if self.eliminate_step(i, j) {
                i += 1;
                if i % 1000 == 0 {
                    let nnz: usize = self.rows.iter().map(|r| r.len()).sum();
                    debug!("  snf progress: {i} pivots ({j}/{n} cols), nnz {nnz}, ops {}", self.ops.len());
                }
            }
        }
    }

    fn eliminate_step(&mut self, i: usize, j: usize) -> bool {
        // pivot: the row with minimal nnz among rows ≥ i having an entry in col j
        let Some(i_p) = self.cols[j].iter()
            .filter(|&&r| r >= i)
            .min_by_key(|&&r| self.rows[r].len())
            .copied()
        else {
            return false
        };

        if i_p > i {
            self.swap_rows(i, i_p);
        }
        if j > i {
            self.swap_cols(i, j);
        }

        let u = self.rows[i][&i].normalizing_unit();
        if !u.is_one() {
            self.mul_col(i, &u);
        }

        self.eliminate_at(i);
        self.rank += 1;

        true
    }

    // Clears row i and column i by Euclidean row/col operations on the pivot (i, i).
    fn eliminate_at(&mut self, i: usize) {
        debug_assert!(self.rows[i].contains_key(&i));

        while self.rows[i].len() > 1 || self.cols[i].len() > 1 {
            let modified = self.eliminate_col(i) | self.eliminate_row(i);
            if !modified {
                panic!("sparse snf: no progress at pivot {i}");
            }
        }
    }

    fn eliminate_col(&mut self, i: usize) -> bool {
        let mut modified = false;
        let targets: Vec<usize> = self.cols[i].iter().copied().filter(|&r| r != i).collect();

        for i1 in targets {
            let (Some(x), Some(y)) = (self.entry(i, i).cloned(), self.entry(i1, i).cloned()) else { continue };
            let (d, s, t) = gcdx(&x, &y);
            let (a, b) = (&x / &d, &y / &d);

            // [ s t][rᵢ ]   [ d ]
            // [-b a][rᵢ₁] = [ 0 ]   (at col i)
            self.row_elem([s, t, -b, a], i, i1);
            modified = true;
        }

        modified
    }

    fn eliminate_row(&mut self, i: usize) -> bool {
        let mut modified = false;
        let targets: Vec<usize> = self.rows[i].keys().copied().filter(|&c| c != i).collect();

        for j1 in targets {
            let (Some(x), Some(y)) = (self.entry(i, i).cloned(), self.entry(i, j1).cloned()) else { continue };
            let (d, s, t) = gcdx(&x, &y);
            let (a, b) = (&x / &d, &y / &d);

            // [cᵢ cⱼ₁][s -b]   [d 0]   (at row i)
            //         [t  a] =
            self.col_elem([s, t, -b, a], i, j1);
            modified = true;
        }

        modified
    }

    // ---- diagonal normalization (divisor chain d₁ | d₂ | …) ----

    fn diag_normalize(&mut self) {
        let r = self.rank;
        if r == 0 {
            return
        }

        'outer: loop {
            for i in 0..r-1 {
                if !self.diag_normalize_step(i) {
                    continue 'outer
                }
            }
            break
        }

        for i in 0..r {
            let u = self.rows[i][&i].normalizing_unit();
            if !u.is_one() {
                self.mul_row(i, &u);
            }
        }
    }

    fn diag_normalize_step(&mut self, i: usize) -> bool {
        let x = self.rows[i][&i].clone();
        let y = self.rows[i + 1][&(i + 1)].clone();

        debug_assert!(!x.is_zero() && !y.is_zero());

        if x.divides(&y) {
            return true
        }

        if y.divides(&x) {
            self.swap_rows(i, i + 1);
            self.swap_cols(i, i + 1);
            return false
        }

        // sx + ty = d, a = x/d, b = y/d:
        // [1   1 ][x   ][s  -b]   [d      ]
        // [-tb sa][   y][t   a] = [   xy/d]
        let (d, s, t) = gcdx(&x, &y);
        let (a, b) = (&x / &d, &y / &d);
        let (tb, sa) = (&t * &b, &s * &a);

        self.row_elem([R::one(), R::one(), -tb, sa], i, i + 1);
        self.col_elem([s, t, -b, a], i, i + 1);

        false
    }

    // ---- primitive operations (matrix + op log) ----

    fn entry(&self, i: usize, j: usize) -> Option<&R> {
        self.rows[i].get(&j)
    }

    fn record(&self) -> bool {
        self.flags.iter().any(|&b| b)
    }

    fn swap_rows(&mut self, i: usize, j: usize) {
        let keys: FxHashSet<usize> = self.rows[i].keys().chain(self.rows[j].keys()).copied().collect();
        self.rows.swap(i, j);
        for k in keys {
            let col = &mut self.cols[k];
            let (a, b) = (col.contains(&i), col.contains(&j));
            if a && !b {
                col.remove(&i);
                col.insert(j);
            } else if b && !a {
                col.remove(&j);
                col.insert(i);
            }
        }
        if self.record() {
            self.ops.push(Op::SwapRows(i, j));
        }
    }

    fn swap_cols(&mut self, i: usize, j: usize) {
        let rows: Vec<usize> = self.cols[i].union(&self.cols[j]).copied().collect();
        for r in rows {
            let vi = self.rows[r].remove(&i);
            let vj = self.rows[r].remove(&j);
            if let Some(v) = vj { self.rows[r].insert(i, v); }
            if let Some(v) = vi { self.rows[r].insert(j, v); }
        }
        self.cols.swap(i, j);
        if self.record() {
            self.ops.push(Op::SwapCols(i, j));
        }
    }

    fn mul_row(&mut self, i: usize, u: &R) {
        debug_assert!(u.is_unit());
        for (_, v) in self.rows[i].iter_mut() {
            *v = &*v * u;
        }
        if self.record() {
            self.ops.push(Op::MulRow(i, u.clone()));
        }
    }

    fn mul_col(&mut self, j: usize, u: &R) {
        debug_assert!(u.is_unit());
        let rows: Vec<usize> = self.cols[j].iter().copied().collect();
        for r in rows {
            let v = self.rows[r].get_mut(&j).unwrap();
            *v = &*v * u;
        }
        if self.record() {
            self.ops.push(Op::MulCol(j, u.clone()));
        }
    }

    // rows (i, j) ← (a·rᵢ + b·rⱼ, c·rᵢ + d·rⱼ), det = 1.
    fn row_elem(&mut self, comps: [R; 4], i: usize, j: usize) {
        let [a, b, c, d] = &comps;
        debug_assert!((a * d - b * c).is_one());

        let ri = std::mem::take(&mut self.rows[i]);
        let rj = std::mem::take(&mut self.rows[j]);
        let keys: FxHashSet<usize> = ri.keys().chain(rj.keys()).copied().collect();

        let (mut ni, mut nj) = (FxHashMap::default(), FxHashMap::default());
        for &k in &keys {
            let x = ri.get(&k);
            let y = rj.get(&k);
            let vi = lin(a, x, b, y);
            let vj = lin(c, x, d, y);

            let col = &mut self.cols[k];
            if vi.is_zero() { col.remove(&i); } else { col.insert(i); ni.insert(k, vi); }
            if vj.is_zero() { col.remove(&j); } else { col.insert(j); nj.insert(k, vj); }
        }
        self.rows[i] = ni;
        self.rows[j] = nj;

        if self.record() {
            self.ops.push(Op::Left(comps, i, j));
        }
    }

    // cols (i, j) ← (a·cᵢ + b·cⱼ, c·cᵢ + d·cⱼ), det = 1.
    fn col_elem(&mut self, comps: [R; 4], i: usize, j: usize) {
        let [a, b, c, d] = &comps;
        debug_assert!((a * d - b * c).is_one());

        let rows: Vec<usize> = self.cols[i].union(&self.cols[j]).copied().collect();
        for r in rows {
            let x = self.rows[r].remove(&i);
            let y = self.rows[r].remove(&j);
            let vi = lin(a, x.as_ref(), b, y.as_ref());
            let vj = lin(c, x.as_ref(), d, y.as_ref());

            if vi.is_zero() { self.cols[i].remove(&r); } else { self.cols[i].insert(r); self.rows[r].insert(i, vi); }
            if vj.is_zero() { self.cols[j].remove(&r); } else { self.cols[j].insert(r); self.rows[r].insert(j, vj); }
        }

        if self.record() {
            self.ops.push(Op::Right(comps, i, j));
        }
    }

    // ---- trans materialization ----

    // Replays the op log over sparse identities:
    //   p, qinv accumulate row-wise; pinv, q accumulate column-wise.
    // Inverse of an elementary [a,b;c,d] (det 1) is [d,-b;-c,a], which as a
    // combine on the opposite side takes comps [d, -c, -b, a].
    fn make_trans(&self) -> [Option<SpMat<R>>; 4] {
        let (m, n) = self.shape;
        let [fp, fpinv, fq, fqinv] = self.flags;

        let mut p    = fp.then(|| VecStore::id(m));
        let mut pinv = fpinv.then(|| VecStore::id(m));
        let mut q    = fq.then(|| VecStore::id(n));
        let mut qinv = fqinv.then(|| VecStore::id(n));

        let inv = |[a, b, c, d]: &[R; 4]| -> [R; 4] {
            [d.clone(), -c, -b, a.clone()]
        };

        for op in &self.ops {
            match op {
                Op::SwapRows(i, j) => {
                    if let Some(p) = &mut p { p.swap(*i, *j); }
                    if let Some(pinv) = &mut pinv { pinv.swap(*i, *j); }
                }
                Op::SwapCols(i, j) => {
                    if let Some(q) = &mut q { q.swap(*i, *j); }
                    if let Some(qinv) = &mut qinv { qinv.swap(*i, *j); }
                }
                Op::MulRow(i, u) => {
                    if let Some(p) = &mut p { p.scale(*i, u); }
                    if let Some(pinv) = &mut pinv { pinv.scale(*i, &u.inv().unwrap()); }
                }
                Op::MulCol(j, u) => {
                    if let Some(q) = &mut q { q.scale(*j, u); }
                    if let Some(qinv) = &mut qinv { qinv.scale(*j, &u.inv().unwrap()); }
                }
                Op::Left(comps, i, j) => {
                    if let Some(p) = &mut p { p.combine(comps, *i, *j); }
                    if let Some(pinv) = &mut pinv { pinv.combine(&inv(comps), *i, *j); }
                }
                Op::Right(comps, i, j) => {
                    if let Some(q) = &mut q { q.combine(comps, *i, *j); }
                    if let Some(qinv) = &mut qinv { qinv.combine(&inv(comps), *i, *j); }
                }
            }
        }

        [
            p.map(|s| s.into_spmat_rows(m)),
            pinv.map(|s| s.into_spmat_cols(m)),
            q.map(|s| s.into_spmat_cols(n)),
            qinv.map(|s| s.into_spmat_rows(n)),
        ]
    }
}

// x/d is kept as the Bezout coefficient shortcut when it is a unit:
// s = (x/d)⁻¹ satisfies s·x = d directly (avoids coefficient growth).
fn gcdx<R>(x: &R, y: &R) -> (R, R, R)
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    let (d, s, t) = EucRing::gcdx(x, y);
    let a = x / &d;
    if a.is_unit() {
        let s = a.inv().unwrap();
        (d, s, R::zero())
    } else {
        (d, s, t)
    }
}

fn lin<R>(a: &R, x: Option<&R>, b: &R, y: Option<&R>) -> R
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    match (x, y) {
        (Some(x), Some(y)) => a * x + b * y,
        (Some(x), None) => a * x,
        (None, Some(y)) => b * y,
        (None, None) => R::zero(),
    }
}

// A sequence of sparse vectors (rows of `p`/`qinv` or columns of `pinv`/`q`)
// supporting the replayed operations.
struct VecStore<R> {
    vecs: Vec<FxHashMap<usize, R>>,
}

impl<R> VecStore<R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    fn id(n: usize) -> Self {
        let vecs = (0..n).map(|i| {
            let mut v = FxHashMap::with_capacity_and_hasher(1, Default::default());
            v.insert(i, R::one());
            v
        }).collect();
        Self { vecs }
    }

    fn swap(&mut self, i: usize, j: usize) {
        self.vecs.swap(i, j);
    }

    fn scale(&mut self, i: usize, u: &R) {
        for (_, v) in self.vecs[i].iter_mut() {
            *v = &*v * u;
        }
    }

    // (vᵢ, vⱼ) ← (a·vᵢ + b·vⱼ, c·vᵢ + d·vⱼ)
    fn combine(&mut self, comps: &[R; 4], i: usize, j: usize) {
        let [a, b, c, d] = comps;
        let vi = std::mem::take(&mut self.vecs[i]);
        let vj = std::mem::take(&mut self.vecs[j]);
        let keys: FxHashSet<usize> = vi.keys().chain(vj.keys()).copied().collect();

        let (mut ni, mut nj) = (FxHashMap::default(), FxHashMap::default());
        for k in keys {
            let x = vi.get(&k);
            let y = vj.get(&k);
            let wi = lin(a, x, b, y);
            let wj = lin(c, x, d, y);
            if !wi.is_zero() { ni.insert(k, wi); }
            if !wj.is_zero() { nj.insert(k, wj); }
        }
        self.vecs[i] = ni;
        self.vecs[j] = nj;
    }

    fn into_spmat_rows(self, n_cols: usize) -> SpMat<R> {
        let m = self.vecs.len();
        SpMat::from_entries((m, n_cols), self.vecs.into_iter().enumerate().flat_map(|(i, row)|
            row.into_iter().map(move |(j, v)| (i, j, v))
        ))
    }

    fn into_spmat_cols(self, n_rows: usize) -> SpMat<R> {
        let n = self.vecs.len();
        SpMat::from_entries((n_rows, n), self.vecs.into_iter().enumerate().flat_map(|(j, col)|
            col.into_iter().map(move |(i, v)| (i, j, v))
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::Zero;
    use yui_core::num::FF2;
    use yui_core::poly::Poly;
    use crate::dense::snf::snf;

    // exercises both routes: 0 forces the sparse elimination, MAX the dense fast-pass.
    fn check_snf<R>(a: &SpMat<R>)
    where R: EucRing, for<'x> &'x R: EucRingOps<R> {
        for max in [0, usize::MAX] {
            let s = sp_snf_with(a, [true; 4], max);
            let res = s.result();

            // diagonal shape
            for (i, j, v) in res.iter_nz() {
                assert!(i == j || v.is_zero(), "non-diagonal entry at ({i}, {j}) (max: {max})");
            }

            // divisor chain
            let fs = s.factors();
            for w in fs.windows(2) {
                assert!(w[0].divides(w[1]), "{} does not divide {} (max: {max})", w[0], w[1]);
            }

            // p * a * q = result, pinv * result * qinv = a
            let (p, pinv, q, qinv) = (s.p().unwrap(), s.pinv().unwrap(), s.q().unwrap(), s.qinv().unwrap());
            assert_eq!(&(&(p * a) * q), res, "p*a*q != result (max: {max})");
            assert_eq!(&(&(pinv * res) * qinv), a, "pinv*result*qinv != a (max: {max})");
        }
    }

    fn check_against_dense<R>(a: &SpMat<R>)
    where R: EucRing, for<'x> &'x R: EucRingOps<R> {
        check_snf(a);

        let d = snf(&a.clone().into_dense(), [false; 4]);
        let df: Vec<&R> = d.factors();
        for max in [0, usize::MAX] {
            let s = sp_snf_with(a, [false; 4], max);
            let sf: Vec<&R> = s.factors();
            assert_eq!(sf, df, "factors differ from dense snf (max: {max})");
        }
    }

    #[test]
    fn snf_int() {
        let a: SpMat<i64> = SpMat::from_row_major((3, 3), [1, 2, 3, 4, 5, 6, 7, 8, 9]);
        check_against_dense(&a);
    }

    #[test]
    fn snf_int_tors() {
        let a: SpMat<i64> = SpMat::from_row_major((5, 5), [
            -20, -7, -27, 2, 29,
            17, 8, 14, -4, -10,
            13, 8, 10, -4, -6,
            -9, -2, -14, 0, 16,
            5, 0, 5, -1, -4
        ]);
        let s = sp_snf(&a, [true; 4]);
        check_snf(&a);
        let fs: Vec<i64> = s.factors().into_iter().cloned().collect();
        assert_eq!(fs, vec![1, 1, 1, 2, 60]);
    }

    #[test]
    fn snf_zero() {
        let a: SpMat<i64> = SpMat::zero((3, 4));
        let s = sp_snf(&a, [true; 4]);
        assert_eq!(s.rank(), 0);
        check_snf(&a);
    }

    #[test]
    fn snf_empty() {
        let a: SpMat<i64> = SpMat::zero((0, 0));
        let s = sp_snf(&a, [true; 4]);
        assert_eq!(s.rank(), 0);
    }

    #[test]
    fn snf_rect() {
        let a: SpMat<i64> = SpMat::from_row_major((2, 4), [2, 4, 6, 8, 3, 5, 7, 9]);
        check_against_dense(&a);
    }

    #[test]
    fn snf_poly_h_multiples() {
        // the KhI-residual shape: all entries divisible by H, no unit pivots
        type P = Poly<'H', FF2>;
        let h = P::variable();
        let o = P::zero();
        let h2 = &h * &h;
        let h3 = &h2 * &h;

        let a: SpMat<P> = SpMat::from_row_major((3, 3), [
            h.clone(), h2.clone(), o.clone(),
            o.clone(), h.clone(), h3.clone(),
            h2.clone(), o.clone(), h.clone(),
        ]);
        check_snf(&a);
    }

    #[test]
    fn snf_rand() {
        let a: SpMat<i64> = SpMat::rand((20, 30), 0.2);
        check_snf(&a);
    }

    #[test]
    fn snf_rand_dense_cmp() {
        for _ in 0..5 {
            let a: SpMat<i64> = SpMat::rand((8, 10), 0.4);
            check_against_dense(&a);
        }
    }
}
