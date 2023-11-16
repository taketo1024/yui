use std::fmt::Display;
use itertools::Itertools;
use num_traits::ToPrimitive;
use super::digits::Digits;

pub fn subscript<I>(i: I) -> String
where I: ToPrimitive {
    let i = i.to_isize().unwrap();

    if i == 0 { 
        return '\u{2080}'.into()
    }

    let (init, i) = if i > 0 { 
        (String::new(), i as usize)
    } else { 
        ('\u{208B}'.into(), -i as usize)
    };

    i.digits().into_iter().fold(init, |mut res, d| {
        let c = char::from_u32( ('\u{2080}' as u32) + (d as u32) ).unwrap();
        res.push(c);
        res
    })
}

pub fn superscript<I>(i: I) -> String 
where I: ToPrimitive {
    let i = i.to_isize().unwrap();

    if i == 0 { 
        return '\u{2070}'.into()
    }

    let (init, i) = if i > 0 { 
        (String::new(), i as usize)
    } else { 
        ('\u{207B}'.into(), -i as usize)
    };

    i.digits().into_iter().fold(init, |mut res, d| {
        let c = match d { 
            1 => '\u{00B9}',
            2 => '\u{00B2}',
            3 => '\u{00B3}',
            _ => char::from_u32(('\u{2070}' as u32) + (d as u32)).unwrap()
        };
        res.push(c);
        res
    })
}

pub fn table<S, I, J, I1, I2, D, F>(head: S, rows: I1, cols: I2, entry: F) -> String
where 
    S: Display,
    I: Display,
    J: Display,
    I1: Iterator<Item = I>,
    I2: Iterator<Item = J>,
    D: Display,
    F: Fn(&I, &J) -> D
{
    use prettytable::*;

    let rows = rows.collect_vec();
    let cols = cols.collect_vec();

    fn row<I>(head: String, cols: I) -> Row
    where I: Iterator<Item = String> { 
        let mut cells = vec![Cell::new(head.as_str())];
        cells.extend(cols.map(|str| Cell::new(str.as_str())));
        Row::new(cells)
    }

    let mut table = Table::new();

    table.set_format(*format::consts::FORMAT_CLEAN);
    table.set_titles(row(
        head.to_string(),
        cols.iter().map(|j| j.to_string() )
    ));

    for i in rows.iter() { 
        table.add_row(row(
            i.to_string(),
            cols.iter().map(|j| format!("{}", entry(i, j)))
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
        let table = table("", 1..=3, 4..=6, |i, j| i * 10 + j);
        let a = "    4   5   6 \n 1  14  15  16 \n 2  24  25  26 \n 3  34  35  36 \n";
        assert_eq!(table, a.to_string());
    }
}
