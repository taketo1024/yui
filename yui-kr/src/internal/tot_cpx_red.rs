use log::info;
use num_traits::Zero;
use yui::{EucRing, EucRingOps};
use yui_homology::utils::ChainReducer;
use yui_homology::{isize2, ChainComplexTrait, DisplayTable, GenericChainComplex2, GridTrait};
use yui_matrix::MatTrait;
use yui_matrix::sparse::SpMat;

use super::tot_cpx::KRTotComplex;

pub struct KRTotComplexReducer<'a, R>
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    complex: &'a KRTotComplex<R>,
    reducer: ChainReducer<isize2, R>,
    pub size_limit: usize
}

impl<'a, R> KRTotComplexReducer<'a, R> 
where R: EucRing, for<'x> &'x R: EucRingOps<R> {
    pub fn new(complex: &'a KRTotComplex<R>) -> Self { 
        let reducer = ChainReducer::new(
            complex.support(), 
            complex.d_deg()
        );
        let size_limit = usize::MAX;
        Self { complex, reducer, size_limit }
    }

    pub fn reduce(&mut self) { 
        info!("reduce C_tot (q: {})..", self.complex.q_deg());

        for idx in self.complex.support() {
            self.setup(idx);
            self.reduce_at(idx, false);
        }

        if self.reducer.is_done() { 
            return;
        }

        info!("second run..");

        for idx in self.complex.support() {
            self.reduce_at(idx, true);
        }
    }

    fn setup(&mut self, idx: isize2) {
        let to_idx = idx + self.complex.d_deg();

        let size = { 
            let m = self.complex.rank(to_idx);
            let n = if let Some(q) = self.reducer.matrix(idx) { 
                q.ncols()
            } else { 
                self.complex.rank(idx)
            };
            (m, n)
        };
            
        if size.0.is_zero() || size.1.is_zero() || usize::max(size.0, size.1) > self.size_limit { 
            let d = SpMat::zero(size);
            self.reducer.set_matrix(idx, d, false);
            return;
        }

        // MEMO: The 'next matrix' will serve as trans-back for the current one. 
        // This is to reduce the input dim without `with_trans = true`.
        
        let d = if let Some(q) = self.reducer.matrix(idx) {
            self.complex.d_matrix_for(idx, &q)
        } else { 
            self.complex.d_matrix(idx)
        };

        self.reducer.set_matrix(idx, d, false);

        if self.complex.is_supported(to_idx) { 
            let m = size.0;
            self.reducer.set_matrix(to_idx, SpMat::id(m), false);
        }
    }

    fn reduce_at(&mut self, idx: isize2, deep: bool) { 
        assert!(self.reducer.matrix(idx).is_some());

        let to_idx = idx + self.complex.d_deg();
        let d = self.reducer.matrix(idx).unwrap();

        if d.is_zero() { 
            return;
        }

        info!("d {idx} -> {to_idx}");
        info!("  size:    {:?}", d.shape());

        if usize::max(d.nrows(), d.ncols()) > self.size_limit {
            info!("  skipped.");
            return;
        }
        
        self.reducer.reduce_at(idx, deep);

        let d = self.reducer.matrix(idx).unwrap();
        info!("  reduced: {:?}", d.shape());
    }

    pub fn into_complex(self) -> GenericChainComplex2<R> { 
        let red = self.reducer.into_complex();

        info!("reduced C (q: {})\n{}", self.complex.q_deg(), red.display_table("h", "v"));
        
        red
    }
}