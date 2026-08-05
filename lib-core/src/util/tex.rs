//! LaTeX rendering helpers, gated on the `tex` feature.
//!
//! See: <https://en.wikipedia.org/wiki/LaTeX>


use itertools::Itertools;
use std::fmt::Display;

/// Types that can be rendered as LaTeX math.
pub trait TeX {
    fn tex_math_symbol() -> String;
    fn tex_string(&self) -> String;
}

/// Render a 2D table as a LaTeX `\begin{tabular}` environment.
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
            let s = s.to_string();
            if s.is_empty() { 
                "$ $".to_string()
            } else {
                format!("${}$", s)
            }
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

    // one cell per column, so a table with no columns emits no `&` and stays valid.
    let row = |head: String, cells: Vec<String>|
        std::iter::once(head).chain(cells).join(" & ") + " \\\\\n";

    let hor = row(
        disp(head, math_mode),
        cols.iter().map(|c| disp(c, math_mode)).collect_vec()
    );

    if hor_at_top { 
        res += &hor;
        res += "\\hline\n";
    }

    for i in rows {
        res += &row(
            disp(&i, math_mode),
            cols.iter().map(|j| disp(entry(&i, j), math_mode)).collect_vec()
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

    #[test]
    fn table_without_columns() {
        // the row must carry one cell per column, or it outruns the `{r|}` spec.
        let entry = |i: &i32, j: &i32| i * 10 + j;
        let none: [i32; 0] = [];

        let t = tex_table("Caption", "i", none, none, entry, true, true);
        assert!(t.contains("{r|}"));
        assert!(!t.contains('&'), "empty table has a stray `&`:\n{t}");

        let t = tex_table("Caption", "i", [1, 2], none, entry, true, true);
        assert!(!t.contains('&'), "column-less rows have a stray `&`:\n{t}");

        // and a normal table is unaffected: head + 2 columns = 2 separators per row.
        let t = tex_table("Caption", "i", [1], [4, 5], entry, true, true);
        assert!(t.lines().all(|l| !l.ends_with("\\\\") || l.matches('&').count() == 2), "{t}");
    }
}