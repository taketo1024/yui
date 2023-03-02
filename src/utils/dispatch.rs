macro_rules! dispatch_ring {
    ($method:ident, $c_type:expr $(, $args:expr)*) => {{
        use num_bigint::BigInt;
        use yui_ratio::Ratio;
        use yui_ff::FF;
        use yui_polynomial::{Poly, MPoly};

        match $c_type { 
            CType::ZPoly      => $method::<Poly<'H', i64>>($($args),*),
            CType::ZPoly_128  => $method::<Poly<'H', i128>>($($args),*),
            CType::ZPoly_big  => $method::<Poly<'H', BigInt>>($($args),*),
            CType::ZMpoly     => $method::<MPoly<'H', i64>>($($args),*),
            CType::ZMpoly_128 => $method::<MPoly<'H', i128>>($($args),*),
            CType::ZMpoly_big => $method::<MPoly<'H', BigInt>>($($args),*),
            CType::QMpoly     => $method::<MPoly<'H', Ratio<i64>>>($($args),*),
            CType::QMpoly_128 => $method::<MPoly<'H', Ratio<i128>>>($($args),*),
            CType::QMpoly_big => $method::<MPoly<'H', Ratio<BigInt>>>($($args),*),
            CType::F2Mpoly    => $method::<MPoly<'H', FF<2>>>($($args),*),
            CType::F3Mpoly    => $method::<MPoly<'H', FF<3>>>($($args),*),
            CType::F5Mpoly    => $method::<MPoly<'H', FF<5>>>($($args),*),
            CType::F7Mpoly    => $method::<MPoly<'H', FF<7>>>($($args),*),
            _ => dispatch_eucring!($method, $c_type, $($args),*)
        }
    }};
}

macro_rules! dispatch_eucring {
    ($method:ident, $c_type:expr $(, $args:expr)*) => {{
        use num_bigint::BigInt;
        use yui_ratio::Ratio;
        use yui_ff::FF;
        use yui_quad_int::{GaussInt, EisenInt};
        use yui_polynomial::Poly;

        match $c_type { 
            CType::Z          => $method::<i64>($($args),*),
            CType::Z_128      => $method::<i128>($($args),*),
            CType::Z_big      => $method::<BigInt>($($args),*),
            CType::Q          => $method::<Ratio<i64>>($($args),*),
            CType::Q_128      => $method::<Ratio<i128>>($($args),*),
            CType::Q_big      => $method::<Ratio<BigInt>>($($args),*),
            CType::F2         => $method::<FF<2>>($($args),*),
            CType::F3         => $method::<FF<3>>($($args),*),
            CType::F5         => $method::<FF<5>>($($args),*),
            CType::F7         => $method::<FF<7>>($($args),*),
            CType::QPoly      => $method::<Poly<'H', Ratio<i64>>>($($args),*),
            CType::QPoly_128  => $method::<Poly<'H', Ratio<i128>>>($($args),*),
            CType::QPoly_big  => $method::<Poly<'H', Ratio<BigInt>>>($($args),*),
            CType::F2Poly     => $method::<Poly<'H', FF<2>>>($($args),*),
            CType::F3Poly     => $method::<Poly<'H', FF<3>>>($($args),*),
            CType::F5Poly     => $method::<Poly<'H', FF<5>>>($($args),*),
            CType::F7Poly     => $method::<Poly<'H', FF<7>>>($($args),*),
            CType::Gauss      => $method::<GaussInt<i64>>($($args),*), 
            CType::Gauss_128  => $method::<GaussInt<i128>>($($args),*), 
            CType::Gauss_big  => $method::<GaussInt<BigInt>>($($args),*), 
            CType::Eisen      => $method::<EisenInt<i64>>($($args),*), 
            CType::Eisen_128  => $method::<EisenInt<i128>>($($args),*), 
            CType::Eisen_big  => $method::<EisenInt<BigInt>>($($args),*), 
            _ => err!("{} is not supported for: {}", stringify!($method), $c_type)
        }
    }};
}

pub(crate) use {dispatch_ring, dispatch_eucring};