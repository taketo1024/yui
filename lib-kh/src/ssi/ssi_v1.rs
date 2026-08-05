//! `ssi` via the bigraded `KhIHomology` over `𝔽₂[H]`: build the homology in degrees `0, 1`,
//! then read the `H`-divisibility of each equivariant Lee class off its coordinate vector.
//! Simple and direct, but it carries a full homology computation with basis tracking.

use itertools::Itertools;
use num_traits::Zero;
use log::info;

use yui_core::num::FF2;
use yui_link::{InvLink, Link};

use crate::khi::KhIHomology;
use crate::util::FastPoly;
use crate::util::calc::div_vec;

type P = FastPoly<'H', FF2>;

pub(super) fn ssi_divisibility_v1(l: &InvLink, reduced: bool) -> (i32, i32) {
    let r = if reduced { 1 } else { 2 };
    let c = P::variable();
    let t = P::zero();

    // bottom..=1: building the cheap low degrees and truncating only at the top is faster than
    // the doubly-truncated `0..=1` slice (which widens to the dense `-1..=2`). Builder clamps the start.
    let kh = KhIHomology::new_partial(l, &c, &t, reduced, Some(-(Link::MAX_CROSSING as isize) ..= 1));

    assert_eq!(kh[0].rank(), r);
    assert_eq!(kh[1].rank(), r);    

    info!("KhI[0]: {}", kh[0]);    
    info!("KhI[1]: {}", kh[1]);    

    let zs = kh.canon_cycles();
    
    assert_eq!(zs.len(), 2 * r);
    for (i, z) in zs.iter().enumerate() {
        let expected = if i < r { 0 } else { 1 };
        assert!(!z.is_zero());
        assert_eq!(z.homogeneous_value(|x| kh.h_deg_of(x)), Some(expected));
    }

    let ds = zs.iter().enumerate().map(|(i, z)| {
        let h = kh.h_deg_of_chain(z);
        let v = kh[h].vectorize_euc(z);
        info!("a[{i}] in Kh[{h}]: ({})", v.clone().into_dense().iter().join(","));
        v
    }).map(|v| 
        div_vec(&v.subvec(0..r), &c).expect("invalid divisibility.")
    ).collect_vec();

    let (d0, d1) = if reduced { 
        (ds[0], ds[1])
    } else { 
        assert_eq!(ds[0], ds[1]);
        assert_eq!(ds[2], ds[3]);
        (ds[0], ds[2])
    };

    assert!(d0 <= d1);

    (d0, d1)
}
