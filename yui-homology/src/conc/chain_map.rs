use std::sync::Arc;

use num_traits::Zero;
use yui::lc::{Gen, Lc};
use yui::{Ring, RingOps};

use crate::{ChainComplexTrait, GridDeg, GridTrait};

use super::ChainComplexBase;

/// Represents a chain map between chain complexes.
pub struct ChainMap<I, X, Y, R>
where 
    I: GridDeg,
    X: Gen, Y: Gen,
    R: Ring, for<'x> &'x R: RingOps<R>
{
    deg: I,
    map: Arc<dyn Fn(I, &Lc<X, R>) -> Lc<Y, R>>,
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

    pub fn check_at(&self, source: &ChainComplexBase<I, X, R>, target: &ChainComplexBase<I, Y, R>, i: I) {
        let d_deg = source.d_deg();
        for x in source.get(i).raw_gens().iter() {
            let x = Lc::from(x.clone());
            let dx = source.d(i, &x);
            let fdx = self.apply(i + d_deg, &dx);
            let fx = self.apply(i, &x);
            let dfx = target.d(i, &fx);
            assert!(dfx == fdx, "df != fd for x = {x}.");
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
}

#[cfg(test)]
mod tests { 
    use crate::GenericChainComplex;

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
}