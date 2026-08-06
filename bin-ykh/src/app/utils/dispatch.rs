//! Macros that pick a concrete coefficient ring from the `-t` / `-c` flags and
//! call the generic implementation with it, reporting the required kind (a ring,
//! a Euclidean ring, or a field) when no match exists.

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
    ($mode:ident, $kind:expr, $app:ident, $method:ident, $args:expr) => {{
        use crate::app::utils::dispatch::*;

        $mode!($app, $method, $args)
        .unwrap_or_else(||
            err!("`-t {} -c {}` does not give {}.", $args.c_type(), $args.c_value, $kind)
        )
    }};
}

macro_rules! dispatch_ring {
    ($app:ident, $method:ident, $args:expr) => {{
        use crate::app::utils::dispatch::*;
        dispatch!(try_ring, "a ring", $app, $method, $args)
    }};
}

macro_rules! dispatch_eucring {
    ($app:ident, $method:ident, $args:expr) => {{
        use crate::app::utils::dispatch::*;
        dispatch!(try_eucring, "a Euclidean ring", $app, $method, $args)
    }};
}

macro_rules! dispatch_field {
    ($app:ident, $method:ident, $args:expr) => {{
        use crate::app::utils::dispatch::*;
        dispatch!(try_field, "a field", $app, $method, $args)
    }};
}

pub(crate) use {dispatch, dispatch_ring, dispatch_eucring, dispatch_field};

// -- internal -- //

macro_rules! try_ring {
    ($app:ident, $method:ident, $args:expr) => {{
        if !$args.is_poly() {
            try_std!($app, $method, $args)
        } else if $args.is_euc_ring() {
            try_euc_poly!($app, $method, $args)
        } else {
            try_noneuc_poly!($app, $method, $args)
        }
    }}
}

macro_rules! try_eucring {
    ($app:ident, $method:ident, $args:expr) => {{
        if !$args.is_poly() {
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

        match $args.c_type() {
            CType::Z     => invoke!(Z,  $app, $method, $args),
            CType::Q     => invoke!(Q,  $app, $method, $args),
            CType::F2    => invoke!(F2, $app, $method, $args),
            CType::F3    => invoke!(F3, $app, $method, $args),
        }
    }}
}

macro_rules! try_field {
    ($app:ident, $method:ident, $args:expr) => {{
        use yui_core::num::{Ratio, FF};

        type Q = Ratio<Int>;
        type F2 = FF<2>;
        type F3 = FF<3>;

        match $args.c_type() {
            CType::Q     => invoke!(Q,  $app, $method, $args),
            CType::F2    => invoke!(F2, $app, $method, $args),
            CType::F3    => invoke!(F3, $app, $method, $args),
            _ => None
        }
    }}
}

macro_rules! try_euc_poly {
    ($app:ident, $method:ident, $args:expr) => {{
        use yui_core::num::Ratio;
        use yui_core::num::FF;
        use yui_kh::util::FastPoly;

        type Q = Ratio<Int>;
        type F2 = FF<2>;
        type F3 = FF<3>;

        let vars = $args.poly_vars();

        match ($args.c_type(), vars) {
            (CType::Q,  PolyVars::H) => invoke!(FastPoly<'H', Q>,  $app, $method, $args),
            (CType::Q,  PolyVars::T) => invoke!(FastPoly<'T', Q>,  $app, $method, $args),
            (CType::F2, PolyVars::H) => invoke!(FastPoly<'H', F2>, $app, $method, $args),
            (CType::F2, PolyVars::T) => invoke!(FastPoly<'T', F2>, $app, $method, $args),
            (CType::F3, PolyVars::H) => invoke!(FastPoly<'H', F3>, $app, $method, $args),
            (CType::F3, PolyVars::T) => invoke!(FastPoly<'T', F3>, $app, $method, $args),
            _ => None
        }
    }}
}

macro_rules! try_noneuc_poly {
    ($app:ident, $method:ident, $args:expr) => {{
        use yui_core::num::Ratio;
        use yui_core::num::FF;
        use yui_core::poly::Poly2;
        use yui_kh::util::FastPoly;

        type Z = Int;
        type Q = Ratio<Int>;
        type F2 = FF<2>;
        type F3 = FF<3>;

        let vars = $args.poly_vars();

        match ($args.c_type(), vars) {
            (CType::Z,  PolyVars::H ) => invoke!(FastPoly<'H', Z>, $app, $method, $args),
            (CType::Z,  PolyVars::T ) => invoke!(FastPoly<'T', Z>, $app, $method, $args),
            (CType::Z,  PolyVars::HT) => invoke!(Poly2<'H', 'T', Z>, $app, $method, $args),
            (CType::Q,  PolyVars::HT) => invoke!(Poly2<'H', 'T', Q>, $app, $method, $args),
            (CType::F2, PolyVars::HT) => invoke!(Poly2<'H', 'T', F2>, $app, $method, $args),
            (CType::F3, PolyVars::HT) => invoke!(Poly2<'H', 'T', F3>, $app, $method, $args),
            _ => None
        }
    }}
}

macro_rules! invoke {
    ($c_type:ty, $app:ident, $method:ident, $args:expr) => {{
        let res = $app::<$c_type>::$method($args);
        Some(res)
    }}
}

pub(crate) use {invoke, try_field, try_ring, try_eucring, try_std, try_euc_poly, try_noneuc_poly};
