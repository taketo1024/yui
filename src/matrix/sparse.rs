use std::ops::Deref;

use num_traits::Zero;
use sprs::{CsMatBase, SpIndex};

pub trait CsMatElem: Clone + Default + Send + Sync + sprs::MulAcc {}

impl<T> CsMatElem for T
    where Self: Clone + Default + Send + Sync + sprs::MulAcc 
{}

pub trait CsMatExt { 
    fn is_zero(&self) -> bool;
}

impl<N, I, IptrStorage, IndStorage, DataStorage, Iptr> CsMatExt for CsMatBase<N, I, IptrStorage, IndStorage, DataStorage, Iptr>
where
    N: Zero,
    I: SpIndex,
    Iptr: SpIndex,
    IptrStorage: Deref<Target = [Iptr]>,
    IndStorage: Deref<Target = [I]>,
    DataStorage: Deref<Target = [N]>,
{
    fn is_zero(&self) -> bool {
        self.data().iter().all(|a| a.is_zero())
    }
}