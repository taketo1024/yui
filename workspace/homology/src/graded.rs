use itertools::Itertools;
use yui_core::{isize2, usize2, Deg};

pub trait Graded<I>
where I: Deg { 
    type Itr: Iterator<Item = I>;

    fn support(&self) -> Self::Itr;
    fn is_supported(&self, i: I) -> bool { 
        self.support().contains(&i)
    }

    fn display(&self, i: I) -> String;
}

pub trait PrintSeq<I> {
    fn display_seq(&self) -> String;
    fn print_seq(&self) {
        println!("{}", self.display_seq())
    }
}

macro_rules! impl_print_seq {
    ($t:ident) => {
        impl<T> PrintSeq<$t> for T
        where T: Graded<$t> {
            fn display_seq(&self) -> String {
                use yui_utils::table;
                let str = table("i", [""].iter(), self.support(), |_, &i| {
                    self.display(i)
                });
                str
            }
        }                
    };
}

impl_print_seq!(isize);
impl_print_seq!(usize);

pub trait PrintTable<I> {
    fn display_table(&self) -> String;
    fn print_table(&self) {
        println!("{}", self.display_table())
    }
}

macro_rules! impl_print_table {
    ($t:ident) => {
        impl<T> PrintTable<$t> for T
        where T: Graded<$t> {
            fn display_table(&self) -> String {
                use yui_utils::table;
        
                let cols = self.support().map(|$t(i, _)| i).unique().sorted();
                let rows = self.support().map(|$t(_, j)| j).unique().sorted().rev();
        
                let str = table("j\\i", rows, cols, |&j, &i| {
                    self.display($t(i, j))
                });
        
                str
            }
        }
    };
}

impl_print_table!(isize2);
impl_print_table!(usize2);

#[cfg(test)]
mod tests { 
    use std::ops::Range;
    use std::vec::IntoIter;
    use super::*;

    struct X(usize);
    impl Graded<usize> for X {
        type Itr = Range<usize>;
        fn support(&self) -> Self::Itr {
            0 .. self.0
        }
        fn display(&self, i: usize) -> String {
            format!("a{i}")
        }
    }

    struct Y(usize, usize);
    impl Graded<usize2> for Y {
        type Itr = IntoIter<usize2>;
        fn support(&self) -> Self::Itr {
            (0 .. self.0).flat_map(|i| 
                (0 .. self.1).map(move |j| usize2(i, j))
            ).collect_vec().into_iter()
        }
        fn display(&self, idx: usize2) -> String {
            let usize2(i, j) = idx;
            format!("a{i}{j}")
        }
    }

    #[test]
    fn seq() { 
        let s = X(3);
        s.print_seq();
    }

    #[test]
    fn table() { 
        let s = Y(3, 2);
        s.print_table();
    }
}