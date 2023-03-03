cfg_if::cfg_if! { 
    if #[cfg(feature = "i128")] { 
        pub type Int = i128;
    } else if #[cfg(feature = "bigint")] { 
        pub type Int = num_bigint::BigInt;
    } else { 
        pub type Int = i64;
    }
}

macro_rules! call {
    ($c_type:ty, $method:ident $(, $args:expr)*) => {
        $method::<$c_type>($($args),*)
    }
}

pub(crate) use call;

macro_rules! dispatch_eucring {
    ($c_type:expr, $method:ident $(, $args:expr)*) => {{
        use crate::utils::dispatch::*;
        use yui_ratio::Ratio;
        use yui_ff::FF;
        use yui_quad_int::{GaussInt, EisenInt};
        use yui_polynomial::Poly;
        type QQ = Ratio<Int>;

        match $c_type { 
            CType::Z      => call!(Int,   $method, $($args),*),
            CType::Q      => call!(QQ,    $method, $($args),*),
            CType::F2     => call!(FF<2>, $method, $($args),*),
            CType::F3     => call!(FF<3>, $method, $($args),*),
            CType::QPoly  => call!(Poly<'H', QQ>,    $method, $($args),*),
            CType::F2Poly => call!(Poly<'H', FF<2>>, $method, $($args),*),
            CType::F3Poly => call!(Poly<'H', FF<3>>, $method, $($args),*),
            CType::Gauss  => call!(GaussInt<Int>, $method, $($args),*), 
            CType::Eisen  => call!(EisenInt<Int>, $method, $($args),*), 
            _ => err!("{} is not supported for: {}", stringify!($method), $c_type)
        }
    }};
}

macro_rules! dispatch_ring {
    ($c_type:expr, $method:ident $(, $args:expr)*) => {{
        use crate::utils::dispatch::*;
        use yui_ratio::Ratio;
        use yui_ff::FF;
        use yui_polynomial::{Poly, MPoly};
        type QQ = Ratio<Int>;

        match $c_type { 
            CType::ZPoly   => call!( Poly<'H', Int>,   $method, $($args),*),
            CType::ZMpoly  => call!(MPoly<'H', Int>,   $method, $($args),*),
            CType::QMpoly  => call!(MPoly<'H', QQ>,    $method, $($args),*),
            CType::F2Mpoly => call!(MPoly<'H', FF<2>>, $method, $($args),*),
            CType::F3Mpoly => call!(MPoly<'H', FF<3>>, $method, $($args),*),
            _ => dispatch_eucring!($c_type, $method, $($args),*)
        }
    }};
}

pub(crate) use {dispatch_ring, dispatch_eucring};