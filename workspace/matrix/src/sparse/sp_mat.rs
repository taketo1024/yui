use std::ops::{Add, AddAssign, Neg, Sub, SubAssign, Mul, MulAssign, Range};
use std::iter::zip;
use std::fmt::Display;
use num_traits::{Zero, One};
use sprs::{TriMat, CsMat, PermView, CsVecView};
use auto_impl_ops::auto_ops;
use yui_core::{Ring, RingOps, AddMonOps, AddGrpOps, MonOps, AddGrp};
use crate::dense::*;
use super::sp_vec::SpVec;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SpMat<R> { 
    cs_mat: CsMat<R>
}

impl<R> MatType for SpMat<R> {
    fn shape(&self) -> (usize, usize) {
        self.cs_mat.shape()
    }
}

impl<R> SpMat<R> { 
    pub fn cs_mat(&self) -> &CsMat<R> { 
        &self.cs_mat
    }

    pub fn cs_mat_into(self) -> CsMat<R> { 
        self.cs_mat
    }

    pub fn iter(&self) -> impl Iterator<Item = (usize, usize, &R)> { 
        self.cs_mat.iter().map(|(a, (i, j))| {
            (i, j, a)
        })
    }
}

impl<R> From<CsMat<R>> for SpMat<R> {
    fn from(cs_mat: CsMat<R>) -> Self {
        assert!(cs_mat.is_csc());
        Self { cs_mat }
    }
}

impl<R> From<&Mat<R>> for SpMat<R>
where R: Clone + Zero {
    fn from(a: &Mat<R>) -> Self {
        let n = a.cols();
        SpMat::generate(a.shape(), |set| { 
            for (k, a) in a.array().iter().enumerate() {
                if a.is_zero() { continue }
                let (i, j) = (k / n, k % n);
                set(i, j, a.clone());
            }
        })
    }
}

impl<R> Display for SpMat<R>
where R: Clone + Zero + Display {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.to_dense().fmt(f)
    }
}

impl<R> Default for SpMat<R>
where R: Clone + Default + Zero {
    fn default() -> Self {
        Self::zero((0, 0))
    }
}

impl<R> SpMat<R>
where R: Clone {
    pub fn col_vec(&self, j: usize) -> SpVec<R> { 
        let cs_vec = self.col_view(j).to_owned();
        SpVec::from(cs_vec)
    }
}

macro_rules! _generate {
    ($shape:expr, $f:expr) => {{
        let mut t = TriMat::new($shape);
        $f( &mut |i, j, a| { 
            if !a.is_zero() { 
                t.add_triplet(i, j, a)
            }
        });
        let cs_mat = t.to_csc();
        Self::from(cs_mat)
    }}
}

impl<R> SpMat<R> where R: Clone + Zero { 
    pub fn generate<F>(shape: (usize, usize), f: F) -> Self
    where F: FnOnce(&mut (dyn FnMut(usize, usize, R))) { 
        _generate!(shape, f)
    }

    pub fn from_vec(shape: (usize, usize), grid: Vec<R>) -> Self { 
        let n = shape.1;
        Self::generate(shape, |init| { 
            grid.into_iter().enumerate().for_each(|(k, a)| { 
                let (i, j) = (k / n, k % n);
                init(i, j, a)
            })
        })
    }

    pub fn from_entries<T: IntoIterator<Item = (usize, usize, R)>>(shape: (usize, usize), iter: T) -> Self {
        Self::generate(shape, |init| { 
            for (i, j, a) in iter { 
                init(i, j, a)
            }
        })
    }

    pub fn transpose(&self) -> Self { 
        self.view().transpose().to_owned()
    }

    pub fn permute<'b>(&self, p: PermView<'b>, q: PermView<'b>) -> SpMat<R> { 
        self.view().permute(p, q).to_owned()
    }

    pub fn permute_rows<'b>(&self, p: PermView<'b>) -> SpMat<R> { 
        let id = PermView::identity(self.cols());
        self.permute(p, id)
    }
    
    pub fn permute_cols<'b>(&self, q: PermView<'b>) -> SpMat<R> { 
        let id = PermView::identity(self.rows());
        self.permute(id, q)
    }

    pub fn submat(&self, rows: Range<usize>, cols: Range<usize>) -> SpMat<R> { 
        self.view().submat(rows, cols).to_owned()
    }

    pub fn submat_rows(&self, rows: Range<usize>) -> SpMat<R> { 
        let n = self.cols();
        self.submat(rows, 0 .. n)
    }

    pub fn submat_cols(&self, cols: Range<usize>) -> SpMat<R> { 
        let m = self.rows();
        self.submat(0 .. m, cols)
    }

    pub fn combine_blocks(blocks: [&SpMat<R>; 4]) -> SpMat<R> {
        let [a, b, c, d] = blocks;

        assert_eq!(a.rows(), b.rows());
        assert_eq!(c.rows(), d.rows());
        assert_eq!(a.cols(), c.cols());
        assert_eq!(b.cols(), d.cols());

        let (m, n) = (a.rows() + c.rows(), a.cols() + b.cols());
        let (k, l) = a.shape();
        let targets = zip(
            [a, b, c, d], 
            [(0,0), (0,l), (k,0), (k,l)]
        );

        Self::generate((m, n), |set| { 
            for (x, (di, dj)) in targets { 
                for (i, j, r) in x.iter() {
                    set(i + di, j + dj, r.clone());
                }
            }    
        })
    }

    pub fn concat(&self, b: &Self) -> Self { 
        let zero = |m, n| SpMat::<R>::zero((m, n));
        Self::combine_blocks([
            self, 
            b, 
            &zero(0, self.cols()), 
            &zero(0, b.cols())
        ])
    }

    pub fn stack(&self, b: &Self) -> Self { 
        let zero = |m, n| SpMat::<R>::zero((m, n));
        Self::combine_blocks([
            self, 
            &zero(self.rows(), 0), 
            b, 
            &zero(b.rows(), 0)
        ])
    }

    pub fn to_dense(&self) -> Mat<R> { 
        self.into()
    }
}

impl<R> SpMat<R>
where R: Clone + Zero + Send + Sync { 
    pub fn generate_sync<F>(shape: (usize, usize), f: F) -> Self
    where F: FnOnce(&mut (dyn FnMut(usize, usize, R) + Send + Sync)) { 
        _generate!(shape, f)
    }
}

impl<R> SpMat<R>
where R: Clone + Zero { 
    pub fn zero(shape: (usize, usize)) -> Self {
        Self::generate(shape, |_| ())
    }

    pub fn is_zero(&self) -> bool {
        self.cs_mat.data().iter().all(|a| a.is_zero())
    }
}

impl<R> SpMat<R>
where R: Clone + Zero + One { 
    pub fn id(n: usize) -> Self { 
        let indptr = (0..=n).collect();
        let indices = (0..n).collect();
        let data = vec![R::one(); n];
        let cs_mat = CsMat::new_csc((n, n), indptr, indices, data);
        Self::from(cs_mat)
    }

    pub fn is_id(&self) -> bool
    where R: PartialEq {
        self.is_square() && self.cs_mat.into_iter().all(|(a, (i, j))| 
            (i == j && a.is_one()) || (i != j && a.is_zero())
        )
    }
}

impl<R> Neg for SpMat<R>
where R: AddGrp, for<'a> &'a R: AddGrpOps<R> {
    type Output = Self;
    fn neg(self) -> Self::Output {
        -&self
    }
}

impl<R> Neg for &SpMat<R>
where R: AddGrp, for<'a> &'a R: AddGrpOps<R> {
    type Output = SpMat<R>;
    fn neg(self) -> Self::Output {
        let neg = self.cs_mat.map(|a| -a);
        SpMat::from(neg)
    }
}

macro_rules! impl_binop {
    ($trait:ident, $method:ident, $r_trait:ident, $r_op_trait:ident) => {
        #[auto_ops]
        impl<'a, 'b, R> $trait<&'b SpMat<R>> for &'a SpMat<R>
        where R: $r_trait, for<'x> &'x R: $r_op_trait<R> {
            type Output = SpMat<R>;
            fn $method(self, rhs: &'b SpMat<R>) -> Self::Output {
                let res = self.cs_mat.$method(&rhs.cs_mat);
                SpMat::from(res)
            }
        }
    };
}

impl_binop!(Add, add, AddGrp, AddGrpOps);
impl_binop!(Sub, sub, AddGrp, AddGrpOps);
impl_binop!(Mul, mul, Ring, RingOps);

macro_rules! impl_ops {
    ($trait:ident, $r_trait:ident, $r_op_trait:ident) => {
        impl<R> $trait<SpMat<R>> for SpMat<R>
        where R: $r_trait, for<'x> &'x R: $r_op_trait<R> {}

        impl<R> $trait<SpMat<R>> for &SpMat<R>
        where R: $r_trait, for<'x> &'x R: $r_op_trait<R> {}
    };
}

impl_ops!(AddMonOps, AddGrp, AddGrpOps);
impl_ops!(AddGrpOps, AddGrp, AddGrpOps);
impl_ops!(MonOps, Ring, RingOps);
impl_ops!(RingOps, Ring, RingOps);

impl<R> SpMat<R> { 
    pub fn view(&self) -> SpMatView<R> { 
        SpMatView::new(self, self.shape(), |i, j| Some((i, j)))
    }

    pub fn col_view(&self, j: usize) -> CsVecView<R> { 
        assert!(j < self.cols());
        self.cs_mat.outer_view(j).unwrap()
    }
}

impl<R> SpMat<R> where R: Clone + Zero + One { 
    // row_perm(p) * a == a.permute_rows(p)
    pub fn from_row_perm(p: PermView) -> Self {
        let n = p.dim();
        Self::generate((n, n), |set| { 
            for i in 0..n { 
                set(p.at(i), i, R::one())
            }
        })
    }

    // a * col_perm(p) == a.permute_cols(p)
    pub fn from_col_perm(p: PermView) -> Self {
        let n = p.dim();
        Self::generate((n, n), |set| { 
            for i in 0..n { 
                set(i, p.at(i), R::one())
            }
        })
    }
}

pub struct SpMatView<'a, 'b, R>  {
    target: &'a SpMat<R>,
    shape: (usize, usize),
    trans: Box<dyn Fn(usize, usize) -> Option<(usize, usize)> + 'b>
}

impl<'a, 'b, R> SpMatView<'a, 'b, R> {
    fn new<F>(target: &'a SpMat<R>, shape: (usize, usize), trans: F) -> Self
    where F: Fn(usize, usize) -> Option<(usize, usize)> + 'b {
        Self { target, shape, trans: Box::new(trans) }
    }

    pub fn shape(&self) -> (usize, usize) { 
        self.shape
    }

    pub fn rows(&self) -> usize { 
        self.shape.0
    }

    pub fn cols(&self) -> usize { 
        self.shape.1
    }

    pub fn iter(&self) -> impl Iterator<Item = (usize, usize, &R)> {
        self.target.iter().filter_map(|(i, j, a)| { 
            (self.trans)(i, j).map(|(i, j)| (i, j, a))
        })
    }

    pub fn transpose(&self) -> SpMatView<R> { 
        SpMatView::new(self.target, (self.shape.1, self.shape.0), move |i, j| (self.trans)(j, i))
    }

    pub fn permute(&self, p: PermView<'b>, q: PermView<'b>) -> SpMatView<R> { 
        SpMatView::new(self.target, self.shape, move |i, j| (self.trans)(p.at(i), q.at(j)))
    }

    pub fn permute_rows(&self, p: PermView<'b>) -> SpMatView<R> { 
        let id = PermView::identity(self.cols());
        self.permute(p, id)
    }
    
    pub fn permute_cols(&self, q: PermView<'b>) -> SpMatView<R> { 
        let id = PermView::identity(self.rows());
        self.permute(id, q)
    }

    pub fn submat(&self, rows: Range<usize>, cols: Range<usize>) -> SpMatView<R> { 
        let (i0, i1) = (rows.start, rows.end);
        let (j0, j1) = (cols.start, cols.end);

        assert!(i0 <= i1 && i1 <= self.rows());
        assert!(j0 <= j1 && j1 <= self.cols());

        SpMatView::new(self.target, (i1 - i0, j1 - j0), move |i, j| { 
            let Some((i, j)) = (self.trans)(i, j) else { 
                return None
            };

            if rows.contains(&i) && cols.contains(&j) {
                Some((i - i0, j - j0))
            } else { 
                None
            }
        })
    }

    pub fn submat_rows(&self, rows: Range<usize>) -> SpMatView<R> { 
        let n = self.cols();
        self.submat(rows, 0 .. n)
    }

    pub fn submat_cols(&self, cols: Range<usize>) -> SpMatView<R> { 
        let m = self.rows();
        self.submat(0 .. m, cols)
    }
}

impl<'a, 'b, R> SpMatView<'a, 'b, R>
where R: Clone + Zero {
    pub fn to_owned(&self) -> SpMat<R> {
        SpMat::generate(self.shape(), |set| { 
            for (i, j, a) in self.iter() { 
                set(i, j, a.clone())
            }
        })
    }

    pub fn to_dense(&self) -> Mat<R> { 
        self.to_owned().to_dense()
    }
}

#[cfg(test)]
impl<R> SpMat<R>
where R: Ring, for<'a> &'a R: RingOps<R> { 
    pub fn rand(shape: (usize, usize), density: f64) -> Self {
        use rand::Rng;
    
        let (m, n) = shape;
        let mut rng = rand::thread_rng();
    
        Self::generate(shape, |set| { 
            for i in 0..m { 
                for j in 0..n { 
                    if rng.gen::<f64>() < density { 
                        set(i, j, R::one());
                    }
                }
            }
        })
    }
}

#[cfg(test)]
pub(super) mod tests { 
    use sprs::PermOwned;

    use super::*;

    #[test]
    fn init() { 
        let a = SpMat::generate((2, 2), |set| { 
            set(0, 0, 1);
            set(0, 1, 2);
            set(1, 0, 3);
            set(1, 1, 4);
        });
        assert_eq!(&a.cs_mat, &CsMat::new_csc((2, 2), vec![0, 2, 4], vec![0, 1, 0, 1], vec![1, 3, 2, 4]));
    }

    #[test]
    fn from_grid() { 
        let a = SpMat::from_vec((2, 2), vec![1,2,3,4]);
        assert_eq!(&a.cs_mat, &CsMat::new_csc((2, 2), vec![0, 2, 4], vec![0, 1, 0, 1], vec![1, 3, 2, 4]));
    }

    #[test]
    fn from_entries() {
        let a = SpMat::from_entries((3,3), vec![(0, 0, 1), (2, 1, 5), (1, 2, 3)]);
        assert_eq!(&a.cs_mat, &CsMat::new_csc((3, 3), vec![0, 1, 2, 3], vec![0, 2, 1], vec![1, 5, 3]));
    }

    #[test]
    fn to_dense() { 
        let a = SpMat::generate((2, 2), |set| { 
            set(0, 0, 1);
            set(0, 1, 2);
            set(1, 0, 3);
            set(1, 1, 4);
        });
        assert_eq!(a.to_dense(), Mat::from(ndarray::array![[1, 2], [3, 4]]));
    }

    #[test]
    fn permute() { 
        let p = PermOwned::new(vec![1,2,3,0]);
        let q = PermOwned::new(vec![3,0,2,1]);
        let a = SpMat::from_vec((4,4), (0..16).collect());
        let b = a.permute(p.view(), q.view());
        assert_eq!(b, SpMat::from_vec((4,4), vec![
            13, 15, 14, 12,
             1,  3,  2,  0,
             5,  7,  6,  4,
             9, 11, 10,  8,
        ]));
    }

    #[test]
    fn transpose() { 
        let a = SpMat::from_vec((3,4), (0..12).collect());
        let b = a.transpose();

        assert_eq!(b, SpMat::from_vec((4,3), vec![
            0, 4, 8, 
            1, 5, 9, 
            2, 6, 10, 
            3, 7, 11, 
        ]));
    }

    #[test]
    fn row_perm() {
        let a = SpMat::from_vec((3, 4), (0..12).collect());
        let p = PermOwned::new(vec![2,0,1]);
        let q = SpMat::from_row_perm(p.view());
        assert!(q * &a == a.permute_rows(p.view()).to_owned())
    }

    #[test]
    fn col_perm() {
        let a = SpMat::from_vec((3, 4), (0..12).collect());
        let p = PermOwned::new(vec![2,0,1,3]);
        let q = SpMat::from_col_perm(p.view());
        assert!(&a * q == a.permute_cols(p.view()).to_owned())
    }
}