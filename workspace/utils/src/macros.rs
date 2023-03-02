#![allow(unused_macros)]
#![allow(unused_imports)]

#[macro_export]
macro_rules! map {
    {$( $key: expr => $val: expr ),*} => {{
        core::convert::From::from([$(($key, $val),)*])
    }}
}

#[macro_export]
macro_rules! set {
    {$($val: expr ),*} => {{
        HashSet::from([$($val,)*])
    }}
}