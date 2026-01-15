use std::sync::OnceLock;

use num_bigint::BigInt;

use crate::dense::Mat;
use crate::dense::lll::{LLLRing, LLLRingOps, lll_hnf};
use crate::utils::dispatcher::{Dispatcher, FnSpec, RefPtr};

static LLL_DISPATCHER: OnceLock<Dispatcher<LLLSpec>> = OnceLock::new();

pub struct LLLSpec;
impl FnSpec for LLLSpec {
    type Input<R: 'static> = (RefPtr<Mat<R>>, [bool; 2]);
    type Output<R: 'static> = (Mat<R>, Option<Mat<R>>, Option<Mat<R>>);
}

pub fn lll_dispatcher() -> &'static Dispatcher<LLLSpec> {
    LLL_DISPATCHER.get_or_init(|| {
        let map = Dispatcher::<LLLSpec>::new();

        fn call_lll_hnf<R>(input: &(RefPtr<Mat<R>>, [bool; 2])) -> (Mat<R>, Option<Mat<R>>, Option<Mat<R>>)
        where R: LLLRing, for<'a> &'a R: LLLRingOps<R> { 
            let (a, flags) = input;
            let a = unsafe { a.as_ref() };
            lll_hnf(a, *flags)
        }

        map.set::<i32>(call_lll_hnf);
        map.set::<i64>(call_lll_hnf);
        map.set::<i128>(call_lll_hnf);
        map.set::<BigInt>(call_lll_hnf);
        
        map
    })
}