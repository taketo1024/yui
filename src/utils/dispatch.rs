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
            err!("method `{}` is not supported for: -c {} -t {}", stringify!($method), $c_value, $c_type)
        )
    }};
}

macro_rules! dispatch_eucring {
    ($c_value:expr, $c_type:expr, $method:ident $(, $args:expr)*) => {{
        use crate::utils::dispatch::*;

        try_eucring!($c_value, $c_type, $method, $($args),*)
        .unwrap_or( 
            err!("method `{}` is not supported for: -c {} -t {}", stringify!($method), $c_value, $c_type)
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

        type Z = Int;
        type Q = Ratio<Int>;
        type F2 = FF<2>;
        type F3 = FF<3>;

        match $c_type {
            CType::Z  => call!(Z,  $method, $($args),*),
            CType::Q  => call!(Q,  $method, $($args),*),
            CType::F2 => call!(F2, $method, $($args),*),
            CType::F3 => call!(F3, $method, $($args),*),
            _         => None
        }
    }}
}

pub(crate) enum PolyVars { 
    H, T, HT
}

pub(crate) fn poly_vars(c_value: &String) -> PolyVars { 
    use std::collections::HashSet;
    
    let s: HashSet<_> = c_value.split(",").collect();
    match (s.contains("H"), s.contains("T")) { 
        (true,  true)  => PolyVars::HT,
        (true,  false) => PolyVars::H,
        (false, true)  => PolyVars::T,
        (false, false) => PolyVars::H
    }
}

macro_rules! try_poly {
    ($c_value:expr, $c_type:expr, $method:ident $(, $args:expr)*) => {{
        cfg_if::cfg_if! {
            if #[cfg(any(feature = "poly", feature = "all"))] {
                use yui_ratio::Ratio;
                use yui_ff::FF;
                use yui_polynomial::Poly;

                type Q = Ratio<Int>;
                type F2 = FF<2>;
                type F3 = FF<3>;

                let vars = poly_vars($c_value);

                match ($c_type, vars) {
                    (CType::QPoly , PolyVars::H) => call!(Poly<'H', Q>,  $method, $($args),*),
                    (CType::QPoly , PolyVars::T) => call!(Poly<'T', Q>,  $method, $($args),*),
                    (CType::F2Poly, PolyVars::H) => call!(Poly<'H', F2>, $method, $($args),*),
                    (CType::F2Poly, PolyVars::T) => call!(Poly<'T', F2>, $method, $($args),*),
                    (CType::F3Poly, PolyVars::H) => call!(Poly<'H', F3>, $method, $($args),*),
                    (CType::F3Poly, PolyVars::T) => call!(Poly<'T', F3>, $method, $($args),*),
                    _ => None
                }
            } else {
                None
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
                use yui_polynomial::{Poly, Poly2};

                type Z = Int;
                type Q = Ratio<Int>;
                type F2 = FF<2>;
                type F3 = FF<3>;

                let vars = poly_vars($c_value);

                match ($c_type, vars) {
                    (CType::ZPoly,  PolyVars::H ) => call!(Poly<'H', Z>, $method, $($args),*),
                    (CType::ZPoly,  PolyVars::T ) => call!(Poly<'T', Z>, $method, $($args),*),
                    (CType::ZPoly,  PolyVars::HT) => call!(Poly2<'H', 'T', Z>, $method, $($args),*),
                    (CType::QPoly,  PolyVars::HT) => call!(Poly2<'H', 'T', Q>, $method, $($args),*),
                    (CType::F2Poly, PolyVars::HT) => call!(Poly2<'H', 'T', F2>, $method, $($args),*),
                    (CType::F3Poly, PolyVars::HT) => call!(Poly2<'H', 'T', F3>, $method, $($args),*),
                    _ => None
                }
            } else {
                None
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
                None
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
