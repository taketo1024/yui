use std::fmt::Display;
use yui_core::{Elem, Ring, RingOps};
use yui_utils::superscript;

pub trait RModStr: Sized + Display
where Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R> { 
    type R;

    fn zero() -> Self;
    fn rank(&self) -> usize;
    fn tors(&self) -> &Vec<Self::R>;

    fn is_free(&self) -> bool { 
        self.tors().is_empty()
    }

    fn is_zero(&self) -> bool { 
        self.rank() == 0 && self.is_free()
    }

    fn fmt_default(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use grouping_by::GroupingBy;
        let rank = self.rank();
        let tors = self.tors();

        if rank == 0 && tors.is_empty() { 
            return f.write_str("0")
        }
    
        let mut res = vec![];
        let symbol = Self::R::set_symbol();
    
        if rank > 1 {
            let str = format!("{}{}", symbol, superscript(rank as isize));
            res.push(str);
        } else if rank == 1 { 
            let str = format!("{}", symbol);
            res.push(str);
        }
    
        for (t, r) in tors.iter().counter(|&t| t.to_string()) { 
            let str = if r > 1 { 
                format!("({}/{}){}", symbol, t, superscript(r as isize))
            } else { 
                format!("({}/{})", symbol, t)
            };
            res.push(str);
        }
    
        let str = res.join(" ⊕ ");
        f.write_str(&str) 
    }
}