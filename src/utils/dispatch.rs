cfg_if::cfg_if! { 
    if #[cfg(feature = "i128")] { 
        pub type Int = i128;
    } else if #[cfg(feature = "bigint")] { 
        pub type Int = num_bigint::BigInt;
    } else { 
        pub type Int = i64;
    }
}

macro_rules! dispatch_ring {
    ($method:ident, $c_type:expr $(, $args:expr)*) => {{
        use crate::utils::dispatch::Int;
        use yui_ratio::Ratio;
        use yui_ff::FF;
        use yui_polynomial::{Poly, MPoly};

        match $c_type { 
            CType::ZPoly      => $method::<Poly<'H', Int>>($($args),*),
            CType::ZMpoly     => $method::<MPoly<'H', Int>>($($args),*),
            CType::QMpoly     => $method::<MPoly<'H', Ratio<Int>>>($($args),*),
            CType::F2Mpoly    => $method::<MPoly<'H', FF<2>>>($($args),*),
            CType::F3Mpoly    => $method::<MPoly<'H', FF<3>>>($($args),*),
            _ => dispatch_eucring!($method, $c_type, $($args),*)
        }
    }};
}

macro_rules! dispatch_eucring {
    ($method:ident, $c_type:expr $(, $args:expr)*) => {{
        use crate::utils::dispatch::Int;
        use yui_ratio::Ratio;
        use yui_ff::FF;
        use yui_quad_int::{GaussInt, EisenInt};
        use yui_polynomial::Poly;

        match $c_type { 
            CType::Z          => $method::<Int>($($args),*),
            CType::Q          => $method::<Ratio<Int>>($($args),*),
            CType::F2         => $method::<FF<2>>($($args),*),
            CType::F3         => $method::<FF<3>>($($args),*),
            CType::QPoly      => $method::<Poly<'H', Ratio<Int>>>($($args),*),
            CType::F2Poly     => $method::<Poly<'H', FF<2>>>($($args),*),
            CType::F3Poly     => $method::<Poly<'H', FF<3>>>($($args),*),
            CType::Gauss      => $method::<GaussInt<Int>>($($args),*), 
            CType::Eisen      => $method::<EisenInt<Int>>($($args),*), 
            _ => err!("{} is not supported for: {}", stringify!($method), $c_type)
        }
    }};
}

pub(crate) use {dispatch_ring, dispatch_eucring};