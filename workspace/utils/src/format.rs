use std::fmt::Display;
use itertools::Itertools;
use super::digits::Digits;

pub fn subscript(i: isize) -> String {
    if i == 0 { return String::from('\u{2080}') }

    let (mut res, i) = if i > 0 { 
        (vec![], i as usize)
    } else { 
        (vec!['\u{208B}'], -i as usize)
    };

    let mut s = i.digits().into_iter().map(|d| {
        char::from_u32(('\u{2080}' as u32) + (d as u32)).unwrap()
    }).collect_vec();

    res.append(&mut s);
    res.iter().collect()
}

pub fn superscript(i: isize) -> String {
    if i == 0 { return String::from('\u{2070}') }

    let (mut res, i) = if i > 0 { 
        (vec![], i as usize)
    } else { 
        (vec!['\u{207B}'], -i as usize)
    };

    let mut s = i.digits().into_iter().map(|d| {
        match d { 
            1 => '\u{00B9}',
            2 => '\u{00B2}',
            3 => '\u{00B3}',
            _ => char::from_u32(('\u{2070}' as u32) + (d as u32)).unwrap()
        }
    }).collect_vec();

    res.append(&mut s);
    res.iter().collect()
}

pub fn table<F, D>(head: &str, rows: &Vec<isize>, cols: &Vec<isize>, entry: F) -> String
where 
    D: Display, 
    F: Fn(isize, isize) -> D
{
    use prettytable::*;

    fn row<I>(head: String, cols: I) -> Row
    where I: Iterator<Item = String> { 
        let mut cells = vec![Cell::new(head.as_str())];
        cells.extend(cols.map(|str| Cell::new(str.as_str())));
        Row::new(cells)
    }

    let mut table = Table::new();

    table.set_format(*format::consts::FORMAT_CLEAN);
    table.set_titles(row(
        String::from(head), 
        cols.iter().map(|j| j.to_string() )
    ));

    for &i in rows { 
        table.add_row(row(
            i.to_string(),
            cols.iter().map(|&j| format!("{}", entry(i, j)))
        ));
    }

    table.to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_subscript() { 
        assert_eq!(subscript(0), "₀");
        assert_eq!(subscript(1234567890), "₁₂₃₄₅₆₇₈₉₀");
        assert_eq!(subscript(-1234567890), "₋₁₂₃₄₅₆₇₈₉₀");
    }

    #[test]
    fn test_superscript() { 
        assert_eq!(superscript(0), "⁰");
        assert_eq!(superscript(1234567890), "¹²³⁴⁵⁶⁷⁸⁹⁰");
        assert_eq!(superscript(-1234567890), "⁻¹²³⁴⁵⁶⁷⁸⁹⁰");
    }

    #[test]
    fn test_table() { 
        let table = table("", &vec![1,2,3], &vec![4,5,6], |i, j| i * 10 + j);
        let a = "    4   5   6 \n 1  14  15  16 \n 2  24  25  26 \n 3  34  35  36 \n";
        assert_eq!(table, a.to_string());
    }
}
