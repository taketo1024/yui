use itertools::Itertools;
use yui_matrix::sparse::*;
use yui_core::{RingOps, Ring};
use crate::{RModStr, RModGrid, GenericChainComplex};

pub trait ChainComplex: RModGrid
where 
    Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R>,
    Self::Output: RModStr<R = Self::R>
{
    fn d_degree(&self) -> Self::Idx;
    fn d_matrix(&self, k: Self::Idx) -> SpMat<Self::R>;

    fn as_generic(&self) -> GenericChainComplex<Self::R, Self::IdxIter> { 
        GenericChainComplex::generate(
            self.indices().clone(), 
            self.d_degree(), 
            |i| Some(self.d_matrix(i))
        )
    }

    fn desc_d_at(&self, i: Self::Idx) -> String {
        let c = |i| self.get(i).map(|v| v.to_string()).unwrap_or("0".to_string());
        let c0 = c(i);
        let c1 = c(i + self.d_degree());
        let d = self.d_matrix(i).to_dense();
        
        format!("C[{i}]: {c0} -> {c1}\n{d}") 
    }

    fn desc_d(&self) -> String { 
        self.indices().filter_map(|i| 
            if !self.is_zero_at(i) || !self.is_zero_at(i + self.d_degree()) {
                Some(self.desc_d_at(i))
            } else { 
                None
            }
        ).join("\n\n")
    }

    fn print_d(&self) {
        println!("{}", self.desc_d());
    }
}