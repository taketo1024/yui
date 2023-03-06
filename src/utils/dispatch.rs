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
    ($c_value:expr, $c_type:expr, $method:ident $(, $args:expr)*) => {{
        use crate::utils::dispatch::*;

        try_ring!($c_value, $c_type, $method, $($args),*)
        .unwrap_or( 
            err!("{} is not supported for: {}", stringify!($method), $c_type) 
        )
    }};
}

macro_rules! dispatch_eucring {
    ($c_value:expr, $c_type:expr, $method:ident $(, $args:expr)*) => {{
        use crate::utils::dispatch::*;

        try_eucring!($c_value, $c_type, $method, $($args),*)
        .unwrap_or( 
            err!("{} is not supported for: {}", stringify!($method), $c_type) 
        )
    }};
}

pub(crate) use {dispatch_ring, dispatch_eucring};

// -- internal -- //

macro_rules! try_ring {
    ($c_value:expr, $c_type:expr, $method:ident $(, $args:expr)*) => {{
        try_eucring!($c_value, $c_type, $method, $($args),*)
        .or_else(|| try_noneuc_poly!($c_value, $c_type, $method, $($args),*) )
    }}
}

macro_rules! try_eucring {
    ($c_value:expr, $c_type:expr, $method:ident $(, $args:expr)*) => {{
        try_std!($c_value, $c_type, $method, $($args),*)
        .or_else(|| try_poly!($c_value, $c_type, $method, $($args),*) )
        .or_else(|| try_qint!($c_value, $c_type, $method, $($args),*) )
    }}
}

macro_rules! try_std {
    ($c_value:expr, $c_type:expr, $method:ident $(, $args:expr)*) => {{
        use yui_ratio::Ratio;
        use yui_ff::FF;

        type Q = Ratio<Int>;
        type F2 = FF<2>;
        type F3 = FF<3>;

        match $c_type {
            CType::Z  => call!(Int, $method, $($args),*),
            CType::Q  => call!(Q,  $method, $($args),*),
            CType::F2 => call!(F2, $method, $($args),*),
            CType::F3 => call!(F3, $method, $($args),*),
            _         => None
        }
    }}
}

macro_rules! try_poly {
    ($c_value:expr, $c_type:expr, $method:ident $(, $args:expr)*) => {{
        cfg_if::cfg_if! {
            if #[cfg(any(feature = "poly", feature = "all"))] {
                use yui_ratio::Ratio;
                use yui_ff::FF;
                use yui_polynomial::Poly;

                type QPoly  = Poly<'H', Ratio<Int>>;
                type F2Poly = Poly<'H', FF<2>>;
                type F3Poly = Poly<'H', FF<3>>;

                match $c_type {
                    CType::QPoly  => call!(QPoly, $method, $($args),*),
                    CType::F2Poly => call!(F2Poly, $method, $($args),*),
                    CType::F3Poly => call!(F3Poly, $method, $($args),*),
                    _             => None
                }
            } else {
                match $c_type {
                    CType::QPoly  |
                    CType::F2Poly |
                    CType::F3Poly => Some(err!("build with `--features poly` to enable polynomial types.")),
                    _             => None
                }
            }
        }
    }}
}

macro_rules! try_noneuc_poly {
    ($c_value:expr, $c_type:expr, $method:ident $(, $args:expr)*) => {{
        cfg_if::cfg_if! {
            if #[cfg(any(feature = "poly", feature = "all"))] {
                use yui_ratio::Ratio;
                use yui_ff::FF;
                use yui_polynomial::{Poly, PolyN};

                type ZPoly  = Poly<'H', Int>;
                type ZMPoly = PolyN<'H', Int>;
                type QMPoly = PolyN<'H', Ratio<Int>>;
                type F2MPoly = PolyN<'H', FF<2>>;
                type F3MPoly = PolyN<'H', FF<3>>;

                match $c_type {
                    CType::ZPoly   => call!(ZPoly,   $method, $($args),*),
                    CType::ZMpoly  => call!(ZMPoly,  $method, $($args),*),
                    CType::QMpoly  => call!(QMPoly,  $method, $($args),*),
                    CType::F2Mpoly => call!(F2MPoly, $method, $($args),*),
                    CType::F3Mpoly => call!(F3MPoly, $method, $($args),*),
                    _              => None
                }
            } else {
                match $c_type {
                    CType::ZPoly   |
                    CType::ZMpoly  |
                    CType::QMpoly  |
                    CType::F2Mpoly |
                    CType::F3Mpoly => Some(err!("build with `--features poly` to enable polynomial types.")),
                    _              => None
                }
            }
        }
    }}
}

macro_rules! try_qint {
    ($c_value:expr, $c_type:expr, $method:ident $(, $args:expr)*) => {{
        cfg_if::cfg_if! {
            if #[cfg(any(feature = "qint", feature = "all"))] {
                type GaussInt = yui_quad_int::GaussInt<Int>;
                type EisenInt = yui_quad_int::GaussInt<Int>;

                match $c_type {
                    CType::Gauss => call!(GaussInt, $method, $($args),*),
                    CType::Eisen => call!(EisenInt, $method, $($args),*),
                    _            => None
                }
            } else {
                match $c_type {
                    CType::Gauss |
                    CType::Eisen => Some(err!("build with `--features qint` to enable quad-int types.")),
                    _            => None
                }
            }
        }
    }}
}

macro_rules! call {
    ($c_type:ty, $method:ident $(, $args:expr)*) => {
        Some($method::<$c_type>($($args),*))
    }
}

pub(crate) use {call, try_ring, try_eucring, try_std, try_poly, try_noneuc_poly, try_qint};
