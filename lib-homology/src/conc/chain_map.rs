use std::sync::Arc;

use num_traits::Zero;
use yui_core::lc::{EitherKey, LcKey, Lc, split_lr};
use yui_core::{EucRing, EucRingOps, Ring, RingOps};
use yui_matrix::sparse::SpMat;

use crate::{GrMod, AddInd, Summand};

use super::ChainComplex;

/// Represents a chain map between chain complexes.
pub struct ChainMap<'a, 'c, I, X, Y, R>
where
    I: AddInd,
    X: LcKey, Y: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    source: &'a ChainComplex<I, X, R>,
    target: &'a ChainComplex<I, Y, R>,
    deg: I,
    map: Arc<dyn Fn(I, &Lc<X, R>) -> Lc<Y, R> + Send + Sync + 'c>,
}

impl<'a, 'c, I, X, Y, R> ChainMap<'a, 'c, I, X, Y, R>
where
    I: AddInd,
    X: LcKey,
    Y: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    pub fn new<F>(
        source: &'a ChainComplex<I, X, R>,
        target: &'a ChainComplex<I, Y, R>,
        deg: I,
        map: F,
    ) -> Self
        where F: Fn(I, &Lc<X, R>) -> Lc<Y, R> + Send + Sync + 'c
    {
        Self { source, target, deg, map: Arc::new(map) }
    }

    pub fn zero(
        source: &'a ChainComplex<I, X, R>,
        target: &'a ChainComplex<I, Y, R>,
        deg: I,
    ) -> Self {
        Self { source, target, deg, map: Arc::new(|_, _| Lc::zero()) }
    }

    pub fn source(&self) -> &'a ChainComplex<I, X, R> {
        self.source
    }

    pub fn target(&self) -> &'a ChainComplex<I, Y, R> {
        self.target
    }

    pub fn deg(&self) -> I {
        self.deg
    }

    pub fn apply(&self, i: I, z: &Lc<X, R>) -> Lc<Y, R> {
        (self.map)(i, z)
    }

    pub fn make_matrix(&self, i: I) -> SpMat<R> {
        self.source[i].make_matrix(&self.target[i + self.deg], |z| self.apply(i, z))
    }

    pub fn make_matrix_euc(&self, i: I) -> SpMat<R>
    where Y: LcKey, R: EucRing, for<'x> &'x R: EucRingOps<R> {
        self.source[i].make_matrix_euc(&self.target[i + self.deg], |z| self.apply(i, z))
    }

    pub fn describe_map(&self) -> String {
        use itertools::Itertools;
        self.source.support().map(|&i| self.describe_map_at(i)).join("\n\n")
    }

    pub fn describe_map_at(&self, i: I) -> String {
        use std::fmt::Write;
        let j = i + self.deg();
        let mut s = format!("({i}) {} -> ({j}) {}", self.source[i], self.target[j]);
        for z in self.source[i].generators() {
            let w = self.apply(i, &z);
            write!(s, "\n\t{z} -> {w}").unwrap();
        }
        s
    }

    pub fn cone<It>(&self, support: It, target_based: bool) -> ChainComplex<I, EitherKey<X, Y>, R>
    where It: IntoIterator<Item = I>, 'c: 'static {
        assert!(self.source.d_deg() == self.target.d_deg());

        let source = self.source;
        let target = self.target;
        let d_deg = source.d_deg();
        let degs = move |i: I| {
            if !target_based {
                (i, i - d_deg)
            } else {
                (i + d_deg, i)
            }
        };

        let summands = GrMod::generate(support, |i| {
            let (i, j) = degs(i);
            let gens = Iterator::chain(
                source[i].raw_generators().iter().map(|x| EitherKey::from_left(x.clone())),
                target[j].raw_generators().iter().map(|y| EitherKey::from_right(y.clone()))
            );
            Summand::from_raw_generators(gens)
        });

        let d1 = source.raw_d();
        let d2 = target.raw_d();
        let f = self.map.clone();

        let d_map = move |i: I, z: &Lc<EitherKey<X, Y>, R>| {
            let (i, j) = degs(i);
            let (x, y) = split_lr(z);

            let dx = d1(i, &x).map_keys(|x2| EitherKey::from_left (x2));
            let fx =  f(i, &x).map_keys(|y2| EitherKey::from_right(y2));
            let dy = d2(j, &y).map_keys(|y2| EitherKey::from_right(y2));

            dx + fx - dy
        };

        ChainComplex::new(summands, d_deg, d_map)
    }

    #[cfg(debug_assertions)]
    pub fn check_for(&self, i: I, x: &X) {
        let d_deg = self.source.d_deg();
        let x = Lc::from(x.clone());
        let dx = self.source.d(i, &x);
        let fdx = self.apply(i + d_deg, &dx);
        let fx = self.apply(i, &x);
        let dfx = self.target.d(i, &fx);

        assert!(dfx == fdx, "df != fd for x = {x}.\n  df = {dfx},\n  fd = {fdx}.");
    }

    #[cfg(debug_assertions)]
    pub fn check_at(&self, i: I) {
        for x in self.source[i].raw_generators().iter() {
            self.check_for(i, x);
        }
    }

    #[cfg(debug_assertions)]
    pub fn check_all(&self) {
        for &i in self.source.support() {
            self.check_at(i);
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{GenericKey, GenericChainComplex1};

    use super::*;

    #[test]
    fn test_s2_to_d3() {
        let s2 = GenericChainComplex1::<i32>::s2();
        let d3 = GenericChainComplex1::<i32>::d3();

        let f = ChainMap::new(&s2, &d3, 0, |_, z| z.clone());

        f.check_all();
    }

    #[test]
    fn test_cone() {
        type T = EitherKey<GenericKey<isize>, GenericKey<isize>>;
        let s2 = GenericChainComplex1::<i32>::s2();
        let d3 = GenericChainComplex1::<i32>::d3();

        let f = ChainMap::new(&s2, &d3, 0, |_, z| z.clone());

        let cone = f.cone((0..=4).rev(), true);
        cone.check_d_all();

        let x = T::from_left(s2[0].raw_generator(0).clone());
        let y = T::from_left(s2[1].raw_generator(0).clone());
        let z = T::from_right(d3[1].raw_generator(0).clone());

        assert_eq!(cone[1].raw_generators().get_index_of(&x), Some(0));
        assert_eq!(cone[2].raw_generators().get_index_of(&y), Some(0));
        assert_eq!(cone[1].raw_generators().get_index_of(&z), Some(4));

        let dx = cone.d(1, &Lc::from(x.clone()));
        assert_eq!(dx, Lc::from(T::from_right(GenericKey(0, 0))));

        let dy = cone.d(2, &Lc::from(y.clone()));
        assert_eq!(dy, Lc::from_iter([
            (T::from_left(GenericKey(0, 0)), -1),
            (T::from_left(GenericKey(0, 1)), 1),
            (T::from_right(GenericKey(1, 0)), 1),
        ]));

        let dz = cone.d(1, &Lc::from(z.clone()));
        assert_eq!(dz, Lc::from_iter([
            (T::from_right(GenericKey(0, 0)), 1),
            (T::from_right(GenericKey(0, 1)), -1),
        ]));
    }
}
