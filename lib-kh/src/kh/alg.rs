use std::fmt::Display;
use std::ops::{Add, AddAssign, Index};
use itertools::join;
use auto_impl_ops::auto_ops;
use num_traits::Zero;
use yui_core::{AddMon, CloneAnd, MathType, Ring, RingOps};
use yui_core::bitseq::{Bit, BitSeq};
use yui_core::lc::{LcKey, Lc};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
pub enum KhAlgGen {
    #[default]
    I,
    X
}

impl KhAlgGen {
    #[allow(non_snake_case)]
    pub fn is_X(&self) -> bool {
        self == &KhAlgGen::X
    }

    pub fn is_1(&self) -> bool {
        self == &KhAlgGen::I
    }

    pub fn deg(&self) -> isize {
        match self {
            KhAlgGen::I => 0,
            KhAlgGen::X => -2
        }
    }

    fn from_bit(b: Bit) -> Self {
        match b {
            Bit::Bit0 => KhAlgGen::I,
            Bit::Bit1 => KhAlgGen::X,
        }
    }

    fn into_bit(self) -> Bit {
        match self {
            KhAlgGen::I => Bit::Bit0,
            KhAlgGen::X => Bit::Bit1
        }
    }
}

impl Display for KhAlgGen {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            KhAlgGen::I => f.write_str("1"),
            KhAlgGen::X => f.write_str("X")
        }
    }
}

impl MathType for KhAlgGen {
    fn math_symbol() -> String {
        "A".to_string()
    }
}

impl LcKey for KhAlgGen {}

#[derive(Clone, Copy, Default, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub struct KhTensor(
    BitSeq
);

impl KhTensor {
    pub fn empty() -> Self {
        Self(BitSeq::empty())
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn iter(&self) -> impl Iterator<Item = KhAlgGen> + use<> {
        self.0.iter().map(KhAlgGen::from_bit
        )
    }

    pub fn set(&mut self, i: usize, x: KhAlgGen) {
        assert!(i < self.len());
        self.0.set(i, x.into_bit());
    }

    pub fn push(&mut self, x: KhAlgGen) {
        self.0.push(x.into_bit());
    }

    pub fn append(&mut self, other: KhTensor) {
        self.0.append(other.0)
    }

    pub fn remove(&mut self, i: usize) {
        self.0.remove(i)
    }

    pub fn insert(&mut self, i: usize, x: KhAlgGen) {
        self.0.insert(i, x.into_bit());
    }

    pub fn sub(&self, l: usize) -> Self {
        Self(self.0.sub(l))
    }

    pub fn is_sub(&self, other: &Self) -> bool {
        self.0.is_sub(&other.0)
    }

    pub fn generate(len: usize) -> impl Iterator<Item = Self> {
        BitSeq::generate(len).map(KhTensor)
    }

    pub fn apply_at<F, R>(&self, i: usize, f: F) -> Lc<KhTensor, R>
    where F: Fn(&KhAlgGen) -> Lc<KhAlgGen, R>, R: Ring, for<'x> &'x R: RingOps<R> {
        assert!(i < self.len());

        f(&self[i]).map_keys(|y|
            self.clone_and(|t| t.set(i, y))
        )
    }

    pub fn apply_each<F, R>(&self, f: F) -> Lc<KhTensor, R>
    where F: Fn(&KhAlgGen) -> Lc<KhAlgGen, R>, R: Ring, for<'x> &'x R: RingOps<R> {
        let l = self.len();
        let init = Lc::from(*self);

        (0..l).fold(init, |res, i| {
            Lc::sum(res.iter().map(|(x, r)|
                x.apply_at(i, &f) * r
            ))
        })
    }
}

impl From<KhAlgGen> for KhTensor {
    fn from(x: KhAlgGen) -> Self {
        Self(BitSeq::from(x.into_bit()))
    }
}

impl<const N: usize> From<[KhAlgGen; N]> for KhTensor {
    fn from(xs: [KhAlgGen; N]) -> Self {
        Self::from_iter(xs)
    }
}

impl FromIterator<KhAlgGen> for KhTensor {
    fn from_iter<I: IntoIterator<Item = KhAlgGen>>(iter: I) -> Self {
        Self(BitSeq::from_iter(iter.into_iter().map(|x|
            x.into_bit()
        )))
    }
}

impl Index<usize> for KhTensor {
    type Output = KhAlgGen;

    fn index(&self, index: usize) -> &Self::Output {
        assert!(index < self.len());
        if self.0[index].is_zero() {
            &KhAlgGen::I
        } else {
            &KhAlgGen::X
        }
    }
}

impl Display for KhTensor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        join(self.iter(), "").fmt(f)
    }
}

impl MathType for KhTensor {
    fn math_symbol() -> String {
        String::from("KhT")
    }
}

impl LcKey for KhTensor {}

#[auto_ops]
impl AddAssign<KhTensor> for KhTensor {
    fn add_assign(&mut self, rhs: Self) {
        self.append(rhs);
    }
}

#[auto_ops]
impl AddAssign<KhAlgGen> for KhTensor {
    fn add_assign(&mut self, x: KhAlgGen) {
        self.push(x);
    }
}

#[derive(Clone)]
pub struct KhAlg<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    h: R,
    t: R
}

impl<R> KhAlg<R>
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn new(h: &R, t: &R) -> Self {
        Self {
            h: h.clone(),
            t: t.clone()
        }
    }

    pub fn h(&self) -> &R {
        &self.h
    }

    pub fn t(&self) -> &R {
        &self.t
    }

    pub fn ht(&self) -> (&R, &R) {
        (&self.h, &self.t)
    }

    pub fn mul(&self, x: KhAlgGen, y: KhAlgGen) -> Lc<KhAlgGen, R> {
        use KhAlgGen::{I, X};
        let (h, t) = (self.h(), self.t());

        match (x, y) {
            (I, I) => Lc::from(I),
            (X, I) | (I, X) => Lc::from(X),
            (X, X) => match (h.is_zero(), t.is_zero()) {
                (true,  true ) => Lc::zero(),
                (false, true ) => Lc::from((X, h.clone())),
                (true,  false) => Lc::from((I, t.clone())),
                (false, false) => Lc::from_iter([(X, h.clone()), (I, t.clone())])
            }
        }
    }

    pub fn mul_tensor(&self, x: &KhTensor, in_index: (usize, usize), out_index: usize) -> Lc<KhTensor, R> {
        assert_ne!(in_index.0, in_index.1);

        let (i, j) = if in_index.0 < in_index.1 {
            in_index
        } else {
            (in_index.1, in_index.0)
        };
        let k = out_index;

        self.mul(x[i], x[j]).map_keys(|a| {
            x.clone_and(|y| {
                y.remove(j);
                y.remove(i);
                y.insert(k, a);
            })
        })
    }

    pub fn comul(&self, x: KhAlgGen) -> Lc<KhTensor, R> {
        use KhAlgGen::{I, X};
        let (h, t) = (self.h(), self.t());
        let tsr = |x, y| KhTensor::from_iter([x, y]);

        match x {
            I => if h.is_zero() {
                Lc::from_iter([
                    (tsr(X, I), R::one()),
                    (tsr(I, X), R::one())
                ])
            } else {
                Lc::from_iter([
                    (tsr(X, I), R::one()),
                    (tsr(I, X), R::one()),
                    (tsr(I, I), -h)
                ])
            },
            X => if t.is_zero() {
                Lc::from(
                    (tsr(X, X), R::one())
                )
            } else {
                Lc::from_iter([
                    (tsr(X, X), R::one()),
                    (tsr(I, I), t.clone())
                ])
            }
        }
    }

    pub fn comul_tensor(&self, x: &KhTensor, in_index: usize, out_index: (usize, usize)) -> Lc<KhTensor, R> {
        assert_ne!(out_index.0, out_index.1);

        let i = in_index;
        let (j, k) = if out_index.0 < out_index.1 {
            out_index
        } else {
            (out_index.1, out_index.0)
        };

        self.comul(x[i]).map_keys(|a|
            x.clone_and(|y| {
                y.remove(i);
                y.insert(j, a[0]);
                y.insert(k, a[1]);
            })
        )
    }

    pub fn sigma(&self, x: &KhAlgGen) -> Lc<KhAlgGen, R> {
        match x {
            KhAlgGen::I => Lc::from((KhAlgGen::I, R::one())),
            KhAlgGen::X => Lc::from_iter([
                (KhAlgGen::X, -R::one()),
                (KhAlgGen::I, self.h().clone())
            ])
        }
    }
}

#[cfg(test)]
pub mod tests {
    use num_traits::Zero;
    use yui_core::lc::Lc;

    use super::{KhAlgGen, KhAlg, KhTensor};

    #[test]
    fn alg_gen() {
        use KhAlgGen::{I, X};
        assert_eq!(I.deg(), 0);
        assert_eq!(X.deg(), -2);
    }

    #[test]
    fn str_prod_kh() {
        use KhAlgGen::{I, X};
        let a = KhAlg::new(&0, &0);
        assert_eq!(a.mul(I, I), Lc::from((I, 1)));
        assert_eq!(a.mul(X, I), Lc::from((X, 1)));
        assert_eq!(a.mul(I, X), Lc::from((X, 1)));
        assert_eq!(a.mul(X, X), Lc::zero());
    }

    #[test]
    fn str_coprod_kh() {
        use KhAlgGen::{I, X};
        let a = KhAlg::new(&0, &0);
        assert_eq!(a.comul(I), Lc::from_iter([
            (KhTensor::from([X, I]), 1),
            (KhTensor::from([I, X]), 1),
        ]));
        assert_eq!(a.comul(X), Lc::from(
            (KhTensor::from_iter([X, X]), 1)
        ));
    }
    #[test]
    fn str_prod_bn() {
        use KhAlgGen::{I, X};
        let a = KhAlg::new(&1, &0);
        assert_eq!(a.mul(I, I), Lc::from((I, 1)));
        assert_eq!(a.mul(X, I), Lc::from((X, 1)));
        assert_eq!(a.mul(I, X), Lc::from((X, 1)));
        assert_eq!(a.mul(X, X), Lc::from((X, 1)));
    }

    #[test]
    fn str_coprod_bn() {
        use KhAlgGen::{I, X};
        let a = KhAlg::new(&1, &0);
        assert_eq!(a.comul(I), Lc::from_iter([
            (KhTensor::from([X, I]), 1),
            (KhTensor::from([I, X]), 1),
            (KhTensor::from([I, I]), -1),
        ]));
        assert_eq!(a.comul(X), Lc::from(
            (KhTensor::from_iter([X, X]), 1)
        ));
    }
    #[test]
    fn str_prod_lee() {
        use KhAlgGen::{I, X};
        let a = KhAlg::new(&0, &1);
        assert_eq!(a.mul(I, I), Lc::from((I, 1)));
        assert_eq!(a.mul(X, I), Lc::from((X, 1)));
        assert_eq!(a.mul(I, X), Lc::from((X, 1)));
        assert_eq!(a.mul(X, X), Lc::from((I, 1)));
    }

    #[test]
    fn str_coprod_lee() {
        use KhAlgGen::{I, X};
        let a = KhAlg::new(&0, &1);
        assert_eq!(a.comul(I), Lc::from_iter([
            (KhTensor::from([X, I]), 1),
            (KhTensor::from([I, X]), 1),
        ]));
        assert_eq!(a.comul(X), Lc::from_iter([
            (KhTensor::from([X, X]), 1),
            (KhTensor::from([I, I]), 1),
        ]));
    }

    #[test]
    fn mul_x_at() {
        use KhAlgGen::{I, X};
        let a = KhAlg::new(&2, &1);
        let x = KhTensor::from_iter([I, X, I]);

        assert_eq!(
            x.apply_at(0, |&x| a.mul(x, X)),
            Lc::from(KhTensor::from_iter([X, X, I]))
        );
        assert_eq!(
            x.apply_at(1, |&x| a.mul(x, X)),
            Lc::from_iter([
                (KhTensor::from([I, X, I]), 2),
                (KhTensor::from([I, I, I]), 1),
            ]
        ));
        assert_eq!(
            x.apply_at(2, |&x| a.mul(x, X)),
            Lc::from(KhTensor::from_iter([I, X, X]))
        );
    }
}
