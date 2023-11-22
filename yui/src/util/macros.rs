#![allow(unused_macros)]
#![allow(unused_imports)]

#[macro_export]
macro_rules! hashmap {
    {$( $key: expr => $val: expr ),*} => {{
        std::collections::HashMap::from_iter([$(($key, $val),)*])
    }}
}

#[macro_export]
macro_rules! set {
    {$($val: expr ),*} => {{
        FromIterator::from_iter([$($val,)*])
    }}
}

pub use {hashmap, set};