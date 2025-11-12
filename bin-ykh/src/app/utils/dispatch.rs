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
    ($app:ident, $args:expr) => {{
        use crate::app::utils::dispatch::*;

        try_ring!($app, $args)
        .unwrap_or_else(|| 
            err!("`{}` is not supported for: -t {} -c {}", stringify!($app), $args.c_type, $args.c_value)
        )
    }};
}

macro_rules! dispatch_eucring {
    ($app:ident, $args:expr) => {{
        use crate::app::utils::dispatch::*;

        try_eucring!($app, $args)
        .unwrap_or( 
            err!("`{}` is not supported for: -t {} -c {}", stringify!($app), $args.c_type, $args.c_value)
        )
    }};
}

pub(crate) use {dispatch_ring, dispatch_eucring};

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
    ($app:ident, $args:expr) => {{
        if poly_vars(&$args.c_value) == PolyVars::None { 
            try_std!($app, $args)
        } else { 
            try_euc_poly!($app, $args)
            .or_else(|| try_noneuc_poly!($app, $args))
        }
    }}
}

macro_rules! try_eucring {
    ($app:ident, $args:expr) => {{
        if poly_vars(&$args.c_value) == PolyVars::None { 
            try_std!($app, $args)
        } else { 
            try_euc_poly!($app, $args)
        }
    }}
}

macro_rules! try_std {
    ($app:ident, $args:expr) => {{
        use yui::num::{Ratio, FF};

        type Z = Int;
        type Q = Ratio<Int>;
        type F2 = FF<2>;
        type F3 = FF<3>;

        match $args.c_type {
            CType::Z     => run!(Z,  $app, $args),
            CType::Q     => run!(Q,  $app, $args),
            CType::F2    => run!(F2, $app, $args),
            CType::F3    => run!(F3, $app, $args),
        }
    }}
}

macro_rules! try_euc_poly {
    ($app:ident, $args:expr) => {{
        cfg_if::cfg_if! {
            if #[cfg(any(feature = "poly", feature = "all"))] {
                use yui::num::Ratio;
                use yui::num::FF;
                use yui::poly::Poly;

                type Q = Ratio<Int>;
                type F2 = FF<2>;
                type F3 = FF<3>;

                let vars = poly_vars(&$args.c_value);

                match ($args.c_type, vars) {
                    (CType::Q,  PolyVars::H) => run!(Poly<'H', Q>,  $app, $args),
                    (CType::Q,  PolyVars::T) => run!(Poly<'T', Q>,  $app, $args),
                    (CType::F2, PolyVars::H) => run!(Poly<'H', F2>, $app, $args),
                    (CType::F2, PolyVars::T) => run!(Poly<'T', F2>, $app, $args),
                    (CType::F3, PolyVars::H) => run!(Poly<'H', F3>, $app, $args),
                    (CType::F3, PolyVars::T) => run!(Poly<'T', F3>, $app, $args),
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
    ($app:ident, $args:expr) => {{
        cfg_if::cfg_if! {
            if #[cfg(any(feature = "poly", feature = "all"))] {
                use yui::num::Ratio;
                use yui::num::FF;
                use yui::poly::{Poly, Poly2};

                type Z = Int;
                type Q = Ratio<Int>;
                type F2 = FF<2>;
                type F3 = FF<3>;

                let vars = poly_vars(&$args.c_value);

                match ($args.c_type, vars) {
                    (CType::Z,  PolyVars::H ) => run!(Poly<'H', Z>, $app, $args),
                    (CType::Z,  PolyVars::T ) => run!(Poly<'T', Z>, $app, $args),
                    (CType::Z,  PolyVars::HT) => run!(Poly2<'H', 'T', Z>, $app, $args),
                    (CType::Q,  PolyVars::HT) => run!(Poly2<'H', 'T', Q>, $app, $args),
                    (CType::F2, PolyVars::HT) => run!(Poly2<'H', 'T', F2>, $app, $args),
                    (CType::F3, PolyVars::HT) => run!(Poly2<'H', 'T', F3>, $app, $args),
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

macro_rules! run {
    ($c_type:ty, $app:ident, $args:expr) => {{
        let mut app: $app<$c_type> = $app::new($args.clone());
        Some(app.run())
    }}
}

pub(crate) use {run, try_ring, try_eucring, try_std, try_euc_poly, try_noneuc_poly};
