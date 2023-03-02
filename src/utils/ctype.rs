use clap::ValueEnum;

#[allow(non_camel_case_types)]
#[derive(Clone, ValueEnum, Debug, derive_more::Display)]
pub enum CType { 
    Z, Z_128, Z_big, 
    Q, Q_128, Q_big, 
    F2, F3, F5, F7,
    ZPoly, ZPoly_128, ZPoly_big, 
    QPoly, QPoly_128, QPoly_big, 
    F2Poly, F3Poly, F5Poly, F7Poly,
    ZMpoly, ZMpoly_128, ZMpoly_big,
    QMpoly, QMpoly_128, QMpoly_big, 
    F2Mpoly, F3Mpoly, F5Mpoly, F7Mpoly,
    Gauss, Gauss_128, Gauss_big,
    Eisen, Eisen_128, Eisen_big,
    None
}