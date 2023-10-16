use std::marker::PhantomData;

use yui_core::{Ring, RingOps};

use super::deg::{Deg, isize2, isize3};

pub type ChainComplex<R>  = ChainComplexBase<isize,  R>;
pub type ChainComplex2<R> = ChainComplexBase<isize2, R>;
pub type ChainComplex3<R> = ChainComplexBase<isize3, R>;

pub struct ChainComplexBase<D, R>
where D: Deg, R: Ring, for<'x> &'x R: RingOps<R> {
    _d: PhantomData<D>,
    _r: PhantomData<R>
}