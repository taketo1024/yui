use std::sync::Arc;

use num_traits::Zero;
use yui_core::lc::{EitherGen, Gen, Lc};
use yui_core::{Ring, RingOps};

use crate::{ChainComplexTrait, Grid, GridDeg, GridTrait, Summand};

use super::ChainComplexBase;

/// Represents a chain map between chain complexes.
pub struct ChainMap<I, X, Y, R>
where 
    I: GridDeg,
    X: Gen, Y: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    deg: I,
    map: Arc<dyn Fn(I, &Lc<X, R>) -> Lc<Y, R> + Send + Sync>,
}

impl<I, X, Y, R> ChainMap<I, X, Y, R>
where 
    I: GridDeg,
    X: Gen,
    Y: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    pub fn new<F>(
        source: &ChainComplexBase<I, X, R>,
        target: &ChainComplexBase<I, Y, R>,
        deg: I,
        map: F,
    ) -> Self
        where F: Fn(I, &Lc<X, R>) -> Lc<Y, R> + Send + Sync + 'static 
    {
        assert!(source.d_deg() == target.d_deg());

        let map = Arc::new(map);
        Self {
            deg,
            map,
        }
    }

    pub fn zero(deg: I
    ) -> Self {
        Self {
            deg,
            map: Arc::new(|_, _| Lc::zero()),
        }
    }

    pub fn deg(&self) -> I {
        self.deg
    }

    pub fn apply(&self, i: I, z: &Lc<X, R>) -> Lc<Y, R> {
        (self.map)(i, &z)
    }

    pub fn check_for(&self, source: &ChainComplexBase<I, X, R>, target: &ChainComplexBase<I, Y, R>, i: I, x: &X) {
        let d_deg = source.d_deg();
        let x = Lc::from(x.clone());
        let dx = source.d(i, &x);
        let fdx = self.apply(i + d_deg, &dx);
        let fx = self.apply(i, &x);
        let dfx = target.d(i, &fx);

        assert!(dfx == fdx, "df != fd for x = {x}.\n  df = {dfx},\n  fd = {fdx}.");
    }


    pub fn check_at(&self, source: &ChainComplexBase<I, X, R>, target: &ChainComplexBase<I, Y, R>, i: I) {
        for x in source.get(i).raw_gens().iter() {
            self.check_for(source, target, i, x);
        }
    }

    pub fn check_all(
        &self,
        source: &ChainComplexBase<I, X, R>,
        target: &ChainComplexBase<I, Y, R>,
    ) {
        for i in source.support() {
            self.check_at(source, target, i);
        }
    }

    pub fn print_map(&self, source: &ChainComplexBase<I, X, R>, target: &ChainComplexBase<I, Y, R>) { 
        for i in source.support() { 
            self.print_map_at(source, target, i);
            println!();
        }
    }

    pub fn print_map_at(&self, source: &ChainComplexBase<I, X, R>, target: &ChainComplexBase<I, Y, R>, i: I) { 
        let j = i + self.deg();
        println!("({i}) {} -> ({j}) {}", source[i], target[j]);
        for z in source[i].gens() { 
            let w = self.apply(i, &z);
            println!("\t{z} -> {w}");
        }
    }

    pub fn cone<It>(&self, source: &ChainComplexBase<I, X, R>, target: &ChainComplexBase<I, Y, R>, support: It, target_based: bool) -> ChainComplexBase<I, EitherGen<X, Y>, R>
    where It: IntoIterator<Item = I> {
        assert!(source.d_deg() == target.d_deg());

        let d_deg = source.d_deg();
        let degs = move |i: I| {
            if !target_based {
                (i, i - d_deg)
            } else {
                (i + d_deg, i)
            }
        };

        let summands = Grid::generate(support, |i| {
            let (i, j) = degs(i);
            let gens = Iterator::chain(
                source.get(i).raw_gens().iter().map(|x| EitherGen::from_left(x.clone())), 
                target.get(j).raw_gens().iter().map(|y| EitherGen::from_right(y.clone()))
            );
            Summand::from_raw_gens(gens)
        });

        let d1 = source.raw_d();
        let d2 = target.raw_d();
        let f = self.map.clone();

        let d_map = move |i: I, z: &Lc<EitherGen<X, Y>, R>| {
            let (i, j) = degs(i);
            let x = z.filter_gens(|x| x.is_left()).map_gens(|x| x.clone().into_left());
            let y = z.filter_gens(|x| x.is_right()).map_gens(|x| x.clone().into_right());

            let dx = d1(i, &x).map_gens(|x2| EitherGen::from_left(x2.clone()));
            let fx = f(i, &x).map_gens(|y2| EitherGen::from_right(y2.clone()));
            let dy = d2(j, &y).map_gens(|y2| EitherGen::from_right(y2.clone()));

            dx + fx - dy
        };

        ChainComplexBase::new(summands, d_deg, d_map)
    }
}

#[cfg(test)]
mod tests { 
    use crate::{EnumGen, GenericChainComplex};

    use super::*;

    #[test]
    fn test_s2_to_d3() { 
        let s2 = GenericChainComplex::<i32>::s2();
        let d3 = GenericChainComplex::<i32>::d3();

        let f = ChainMap::new(
            &s2, 
            &d3, 
            0, 
            |_, z| z.clone()
        );

        f.check_all(&s2, &d3);
    }

    #[test]
    fn test_cone() { 
        type T = EitherGen<EnumGen<isize>, EnumGen<isize>>;
        let s2 = GenericChainComplex::<i32>::s2();
        let d3 = GenericChainComplex::<i32>::d3();

        let f = ChainMap::new(
            &s2, 
            &d3, 
            0, 
            |_, z| z.clone()
        );

        let cone = f.cone(&s2, &d3, (0..=4).rev(), true);
        cone.check_d_all();

        let x = T::from_left(s2[0].raw_gen(0).clone());
        let y = T::from_left(s2[1].raw_gen(0).clone());
        let z = T::from_right(d3[1].raw_gen(0).clone());

        assert_eq!(cone[1].raw_gens().index_of(&x), Some(0));
        assert_eq!(cone[2].raw_gens().index_of(&y), Some(0));
        assert_eq!(cone[1].raw_gens().index_of(&z), Some(4));

        let dx = cone.d(1, &Lc::from(x.clone()));
        assert_eq!(dx, Lc::from(T::from_right(EnumGen(0, 0))));

        let dy = cone.d(2, &Lc::from(y.clone()));
        assert_eq!(dy, Lc::from_iter([
            (T::from_left(EnumGen(0, 0)), -1),
            (T::from_left(EnumGen(0, 1)), 1),
            (T::from_right(EnumGen(1, 0)), 1),
        ]));

        let dz = cone.d(1, &Lc::from(z.clone()));
        assert_eq!(dz, Lc::from_iter([
            (T::from_right(EnumGen(0, 0)), 1),
            (T::from_right(EnumGen(0, 1)), -1),
        ]));
    }
}