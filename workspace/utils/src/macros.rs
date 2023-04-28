#![allow(unused_macros)]
#![allow(unused_imports)]

#[macro_export]
macro_rules! map {
    {$( $key: expr => $val: expr ),*} => {{
        FromIterator::from_iter([$(($key, $val),)*])
    }}
}

#[macro_export]
macro_rules! set {
    {$($val: expr ),*} => {{
        FromIterator::from_iter([$($val,)*])
    }}
}