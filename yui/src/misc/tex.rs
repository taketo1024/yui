#![cfg(feature = "tex")]

use itertools::Itertools;
use std::fmt::Display;
pub trait TeX { 
    fn tex_math_symbol() -> String;
    fn tex_string(&self) -> String;
}

pub fn tex_table<S, I, J, I1, I2, D, F>(caption: &str, head: S, rows: I1, cols: I2, entry: F, math_mode: bool, hor_at_top: bool) -> String
where 
    S: Display,
    I: Display,
    J: Display,
    I1: IntoIterator<Item = I>,
    I2: IntoIterator<Item = J>,
    D: Display,
    F: Fn(&I, &J) -> D
{
    fn disp<S>(s: S, math_mode: bool) -> String where S: Display { 
        if math_mode {
            format!("${}$", s.to_string())
        } else { 
            s.to_string()
        }
    }

    let cols = cols.into_iter().collect_vec();
    let mut res = String::new();
    
    res += r#"\begin{table}
\centering
\begin{tabular}"#;

    res += &format!("{{r|{}}}\n", "l".repeat(cols.len()));

    let hor = format!("{} & {} \\\\\n", 
        disp(head, math_mode), 
        cols.iter().map(|c| disp(c, math_mode)).join(" & ")
    );

    if hor_at_top { 
        res += &hor;
        res += "\\hline\n";
    }

    for i in rows { 
        res += &format!("{} & {} \\\\\n", 
            disp(&i, math_mode), 
            cols.iter().map(|j| disp(entry(&i, j), math_mode)).join(" & ")
        );
    }

    if !hor_at_top { 
        res += "\\hline\n";
        res += &hor;
    }

    res += "\\end{tabular}\n";
    res += &format!("\\caption{{{caption}}}\n");
    res += "\\end{table}\n";
    res
}

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn test_tex_table() { 
        let _table = tex_table("Caption", "i, j", [1, 2, 3], [4, 5, 6, 7], |i, j| i * 10 + j, true, false);
        // println!("{_table}");
    }
}