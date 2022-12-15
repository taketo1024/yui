#![allow(unused_macros)]
#![allow(unused_imports)]

macro_rules! hashmap {
    {$( $key: expr => $val: expr ),*} => {{
        HashMap::from([$(($key, $val),)*])
    }}
}

macro_rules! hashset {
    {$($val: expr ),*} => {{
        HashSet::from([$($val,)*])
    }}
}

pub(crate) use hashmap;
pub(crate) use hashset;