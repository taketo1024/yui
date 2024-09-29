#![allow(unused_macros)]
#![allow(unused_imports)]

#[macro_export]
macro_rules! hashmap {
    {$( $key: expr => $val: expr ),*} => {{
        std::collections::HashMap::<_, _, std::hash::RandomState>::from_iter([$(($key, $val),)*])
    }}
}

#[macro_export]
macro_rules! hashset {
    {$($val: expr ),*} => {{
        std::collections::HashSet::<_, std::hash::RandomState>from_iter([$($val,)*])
    }}
}

pub use {hashmap, hashset};