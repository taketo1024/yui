use std::collections::HashSet;
use std::ops::Index;
use itertools::Itertools;

use yui_core::{Ring, RingOps};
use crate::{RModStr, RModGrid, Idx2};

pub trait PrintTable {
    fn table(&self) -> String;
    fn print_table(&self);
}

impl<T> PrintTable for T
where
    T: Index<Idx2>,
    T: RModGrid<Idx = Idx2>,
    T::R: Ring, for<'x> &'x T::R: RingOps<T::R>,
    <T as Index<Idx2>>::Output: RModStr<R = T::R>
{
    fn table(&self) -> String {
        use yui_utils::table as f_table;

        fn collect<R, F>(range: R, f: F) -> Vec<isize>
        where R: Iterator<Item = Idx2>, F: Fn(Idx2) -> isize { 
            range.map(f).collect::<HashSet<_>>().into_iter().sorted().collect_vec()
        }
        let is = collect(self.indices(), |idx| idx.0);
        let js = collect(self.indices(), |idx| idx.1).into_iter().rev().collect();

        f_table("j\\i", &js, &is, |j, i| {
            let s = &self[Idx2(i, j)];
            if !s.is_zero() {
                format!("{}", s)
            } else { 
                String::from("")
            }
        })
    }

    fn print_table(&self) {
        println!("{}", self.table())
    }
}