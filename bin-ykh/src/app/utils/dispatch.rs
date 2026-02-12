cfg_if::cfg_if! {
    if #[cfg(feature = "i128")] {
        pub type Int = i128;
    } else if #[cfg(feature = "bigint")] {
        pub type Int = num_bigint::BigInt;
    } else {
        pub type Int = i64;
    }
}

macro_rules! dispatch {
    ($mode:ident, $app:ident, $method:ident, $args:expr) => {{
        use crate::app::utils::dispatch::*;

        $mode!($app, $method, $args)
        .unwrap_or_else(|| 
            err!("`{}::{}` is not supported for: -t {} -c {}", stringify!($app), stringify!($method), $args.c_type, $args.c_value)
        )
    }};
}

macro_rules! dispatch_ring {
    ($app:ident, $method:ident, $args:expr) => {{
        use crate::app::utils::dispatch::*;
        dispatch!(try_ring, $app, $method, $args)
    }};
}

macro_rules! dispatch_eucring {
    ($app:ident, $method:ident, $args:expr) => {{
        use crate::app::utils::dispatch::*;
        dispatch!(try_eucring, $app, $method, $args)
    }};
}

pub(crate) use {dispatch, dispatch_ring, dispatch_eucring};

// -- internal -- //

#[derive(PartialEq, Eq)]
pub(crate) enum PolyVars { 
    H, T, HT, None
}

pub(crate) fn poly_vars(c_value: &String) -> PolyVars { 
    use std::collections::HashSet;
    
    let s: HashSet<_> = c_value.split(',').collect();
    match (s.contains("H"), s.contains("T")) { 
        (true,  true)  => PolyVars::HT,
        (true,  false) => PolyVars::H,
        (false, true)  => PolyVars::T,
        (false, false) => PolyVars::None
    }
}

macro_rules! try_ring {
    ($app:ident, $method:ident, $args:expr) => {{
        if poly_vars(&$args.c_value) == PolyVars::None { 
            try_std!($app, $method, $args)
        } else { 
            try_euc_poly!($app, $method, $args)
            .or_else(|| try_noneuc_poly!($app, $method, $args))
        }
    }}
}

macro_rules! try_eucring {
    ($app:ident, $method:ident, $args:expr) => {{
        if poly_vars(&$args.c_value) == PolyVars::None { 
            try_std!($app, $method, $args)
        } else { 
            try_euc_poly!($app, $method, $args)
        }
    }}
}

macro_rules! try_std {
    ($app:ident, $method:ident, $args:expr) => {{
        use yui_core::num::{Ratio, FF};

        type Z = Int;
        type Q = Ratio<Int>;
        type F2 = FF<2>;
        type F3 = FF<3>;

        match $args.c_type {
            CType::Z     => invoke!(Z,  $app, $method, $args),
            CType::Q     => invoke!(Q,  $app, $method, $args),
            CType::F2    => invoke!(F2, $app, $method, $args),
            CType::F3    => invoke!(F3, $app, $method, $args),
        }
    }}
}

macro_rules! try_euc_poly {
    ($app:ident, $method:ident, $args:expr) => {{
        cfg_if::cfg_if! {
            if #[cfg(any(feature = "poly", feature = "all"))] {
                use yui_core::num::Ratio;
                use yui_core::num::FF;
                use yui_core::poly::Poly;

                type Q = Ratio<Int>;
                type F2 = FF<2>;
                type F3 = FF<3>;

                let vars = poly_vars(&$args.c_value);

                match ($args.c_type, vars) {
                    (CType::Q,  PolyVars::H) => invoke!(Poly<'H', Q>,  $app, $method, $args),
                    (CType::Q,  PolyVars::T) => invoke!(Poly<'T', Q>,  $app, $method, $args),
                    (CType::F2, PolyVars::H) => invoke!(Poly<'H', F2>, $app, $method, $args),
                    (CType::F2, PolyVars::T) => invoke!(Poly<'T', F2>, $app, $method, $args),
                    (CType::F3, PolyVars::H) => invoke!(Poly<'H', F3>, $app, $method, $args),
                    (CType::F3, PolyVars::T) => invoke!(Poly<'T', F3>, $app, $method, $args),
                    _ => None
                }
            } else {
                match $c_type {
                    CType::Q  |
                    CType::F2 |
                    CType::F3 => Some(err!("build with `--features poly` to enable polynomial types.")),
                    _         => None
                }
            }
        }
    }}
}

macro_rules! try_noneuc_poly {
    ($app:ident, $method:ident, $args:expr) => {{
        cfg_if::cfg_if! {
            if #[cfg(any(feature = "poly", feature = "all"))] {
                use yui_core::num::Ratio;
                use yui_core::num::FF;
                use yui_core::poly::{Poly, Poly2};

                type Z = Int;
                type Q = Ratio<Int>;
                type F2 = FF<2>;
                type F3 = FF<3>;

                let vars = poly_vars(&$args.c_value);

                match ($args.c_type, vars) {
                    (CType::Z,  PolyVars::H ) => invoke!(Poly<'H', Z>, $app, $method, $args),
                    (CType::Z,  PolyVars::T ) => invoke!(Poly<'T', Z>, $app, $method, $args),
                    (CType::Z,  PolyVars::HT) => invoke!(Poly2<'H', 'T', Z>, $app, $method, $args),
                    (CType::Q,  PolyVars::HT) => invoke!(Poly2<'H', 'T', Q>, $app, $method, $args),
                    (CType::F2, PolyVars::HT) => invoke!(Poly2<'H', 'T', F2>, $app, $method, $args),
                    (CType::F3, PolyVars::HT) => invoke!(Poly2<'H', 'T', F3>, $app, $method, $args),
                    _ => None
                }
            } else {
                match $c_type {
                    CType::Z  |
                    CType::Q  |
                    CType::F2 |
                    CType::F3 => Some(err!("build with `--features poly` to enable polynomial types.")),
                    _         => None
                }
            }
        }
    }}
}

macro_rules! invoke {
    ($c_type:ty, $app:ident, $method:ident, $args:expr) => {{
        let res = $app::<$c_type>::$method($args);
        Some(res)
    }}
}

pub(crate) use {invoke, try_ring, try_eucring, try_std, try_euc_poly, try_noneuc_poly};
