use std::ops::{Add, AddAssign, Neg, Sub, SubAssign, Mul, MulAssign, Range};
use std::iter::zip;
use std::fmt::Display;
use sprs::{TriMat, CsMat, PermView, CsVecView};
use auto_impl_ops::auto_ops;
use crate::math::traits::{Ring, RingOps, AddMonOps, AddGrpOps, MonOps};
use super::DnsMat;
use super::sp_vec::SpVec;

pub use super::dense::MatType;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SpMat<R>
where R: Ring, for<'a> &'a R: RingOps<R> { 
    cs_mat: CsMat<R>
}

impl<R> SpMat<R>
where R: Ring, for<'a> &'a R: RingOps<R> { 
    pub fn generate<F>(shape: (usize, usize), f: F) -> Self
    where F: FnOnce(&mut (dyn FnMut(usize, usize, R) + Send + Sync)) { 
        let mut t = TriMat::new(shape);
        f( &mut |i, j, a| { 
            if !a.is_zero() { 
                t.add_triplet(i, j, a)
            }
        });
        let cs_mat = t.to_csc();
        Self::from(cs_mat)
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

    pub fn to_dense(&self) -> DnsMat<R> { 
        self.into()
    }

    pub fn col_view(&self, j: usize) -> CsVecView<R> { 
        assert!(j < self.cols());
        self.cs_mat.outer_view(j).unwrap()
    }

    pub fn col_vec(&self, j: usize) -> SpVec<R> { 
        let cs_vec = self.col_view(j).to_owned();
        SpVec::from(cs_vec)
    }
}

impl<R> From<CsMat<R>> for SpMat<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    fn from(cs_mat: CsMat<R>) -> Self {
        assert!(cs_mat.is_csc());
        Self { cs_mat }
    }
}

impl<R> From<&DnsMat<R>> for SpMat<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    fn from(a: &DnsMat<R>) -> Self {
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
where R: Ring, for<'a> &'a R: RingOps<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.to_dense().fmt(f)
    }
}

impl<R> Default for SpMat<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    fn default() -> Self {
        Self::zero((0, 0))
    }
}

impl<R> Neg for SpMat<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    type Output = Self;
    fn neg(self) -> Self::Output {
        -&self
    }
}

impl<R> Neg for &SpMat<R>
where R: Ring, for<'a> &'a R: RingOps<R> {
    type Output = SpMat<R>;
    fn neg(self) -> Self::Output {
        let neg = self.cs_mat.map(|a| -a);
        SpMat::from(neg)
    }
}

macro_rules! impl_op {
    ($trait:ident, $method:ident) => {
        #[auto_ops]
        impl<'a, 'b, R> $trait<&'b SpMat<R>> for &'a SpMat<R>
        where R: Ring, for<'x> &'x R: RingOps<R> {
            type Output = SpMat<R>;
            fn $method(self, rhs: &'b SpMat<R>) -> Self::Output {
                let res = self.cs_mat.$method(&rhs.cs_mat);
                SpMat::from(res)
            }
        }
    };
}

impl_op!(Add, add);
impl_op!(Sub, sub);
impl_op!(Mul, mul);

macro_rules! impl_ops {
    ($trait:ident) => {
        impl<R> $trait<SpMat<R>> for SpMat<R>
        where R: Ring, for<'x> &'x R: RingOps<R> {}

        impl<R> $trait<SpMat<R>> for &SpMat<R>
        where R: Ring, for<'x> &'x R: RingOps<R> {}
    };
}

impl_ops!(AddMonOps);
impl_ops!(AddGrpOps);
impl_ops!(MonOps);
impl_ops!(RingOps);

impl<R> MatType for SpMat<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    type R = R;

    fn shape(&self) -> (usize, usize) {
        self.cs_mat.shape()
    }

    fn zero(shape: (usize, usize)) -> Self {
        Self::from(CsMat::zero(shape).to_csc())
    }

    fn is_zero(&self) -> bool {
        self.cs_mat.data().iter().all(|a| a.is_zero())
    }

    fn id(n: usize) -> Self { 
        let indptr = (0..=n).collect();
        let indices = (0..n).collect();
        let data = vec![R::one(); n];
        let cs_mat = CsMat::new_csc((n, n), indptr, indices, data);
        Self::from(cs_mat)
    }

    fn is_id(&self) -> bool {
        self.is_square() && self.cs_mat.into_iter().all(|(a, (i, j))| 
            (i == j && a.is_one()) || (i != j && a.is_zero())
        )
    }
}

impl<R> SpMat<R>
where R: Ring, for<'a> &'a R: RingOps<R> { 
    pub fn view(&self) -> SpMatView<R> { 
        SpMatView::new(self, self.shape(), |i, j| Some((i, j)))
    }

    pub fn transpose(&self) -> SpMatView<R> { 
        SpMatView::new(self, (self.cols(), self.rows()), |i, j| Some((j, i)))
    }

    pub fn permute<'b>(&self, p: PermView<'b>, q: PermView<'b>) -> SpMatView<'_, 'b, R> { 
        SpMatView::new(self, self.shape(), move |i, j| Some((p.at(i), q.at(j))))
    }

    pub fn permute_rows<'b>(&self, p: PermView<'b>) -> SpMatView<'_, 'b, R> { 
        let id = PermView::identity(self.cols());
        self.permute(p, id)
    }
    
    pub fn permute_cols<'b>(&self, q: PermView<'b>) -> SpMatView<'_, 'b, R> { 
        let id = PermView::identity(self.rows());
        self.permute(id, q)
    }

    pub fn submat(&self, rows: Range<usize>, cols: Range<usize>) -> SpMatView<R> { 
        let (i0, i1) = (rows.start, rows.end);
        let (j0, j1) = (cols.start, cols.end);

        assert!(i0 <= i1 && i1 <= self.rows());
        assert!(j0 <= j1 && j1 <= self.cols());

        SpMatView::new(self, (i1 - i0, j1 - j0), move |i, j| { 
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
}

pub struct SpMatView<'a, 'b, R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    target: &'a SpMat<R>,
    shape: (usize, usize),
    trans: Box<dyn Fn(usize, usize) -> Option<(usize, usize)> + 'b>
}

impl<'a, 'b, R> SpMatView<'a, 'b, R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
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

    pub fn to_owned(&self) -> SpMat<R> {
        SpMat::generate(self.shape(), |set| { 
            for (i, j, a) in self.iter() { 
                set(i, j, a.clone())
            }
        })
    }

    pub fn to_dense(&self) -> DnsMat<R> { 
        self.to_owned().to_dense()
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

#[cfg(test)]
pub(super) mod tests { 
    use super::*;

    pub fn mat_rand<R>(shape: (usize, usize), density: f64) -> SpMat<R>
    where R: Ring, for<'a> &'a R: RingOps<R> { 
        use rand::Rng;

        let (m, n) = shape;
        let mut rng = rand::thread_rng();

        SpMat::generate(shape, |set| { 
            for i in 0..m { 
                for j in 0..n { 
                    if rng.gen::<f64>() < density { 
                        set(i, j, R::one());
                    }
                }
            }
        })
    }

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
    fn to_dense() { 
        let a = SpMat::generate((2, 2), |set| { 
            set(0, 0, 1);
            set(0, 1, 2);
            set(1, 0, 3);
            set(1, 1, 4);
        });
        assert_eq!(a.to_dense(), DnsMat::from(ndarray::array![[1, 2], [3, 4]]));
    }
}