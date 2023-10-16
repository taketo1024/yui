use std::fmt::Display;
use itertools::Itertools;

use crate::{Idx2, Grid};

pub trait FmtList { 
    fn list_string(&self) -> String;
    
    fn print_list(&self) { 
        println!("{}", self.list_string())
    }

    fn fmt_list(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result { 
        self.list_string().fmt(f)
    }
}

impl<T> FmtList for T
where
    T: Grid,
    T::Output: Display
{
    fn list_string(&self) -> String {
        let mut res = String::new();
        for (i, v) in self.iter() { 
            res += format!("[{}]: {}\n", i, v).as_str();
        }
        res
    }
}

pub trait DisplayForTable: Display { 
    fn should_display(&self) -> bool { true }
}

pub trait FmtTable {
    fn table_string(&self) -> String;
    
    fn print_table(&self) { 
        println!("{}", self.table_string())
    }

    fn fmt_table(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result { 
        self.table_string().fmt(f)
    }
}

impl<T> FmtTable for T
where
    T: Grid<Idx = Idx2>,
    T::Output: DisplayForTable
{
    fn table_string(&self) -> String {
        use yui_utils::table as f_table;

        let is = self.indices().map(|idx| idx.0).unique().sorted();
        let js = self.indices().map(|idx| idx.1).unique().sorted().rev();

        f_table("j\\i", js, is, |&j, &i| {
            let Some(s) = self.get(Idx2(i, j)) else { 
                return String::from("")
            };
            if s.should_display() {
                format!("{}", s)
            } else { 
                String::from("")
            }
        })
    }
}