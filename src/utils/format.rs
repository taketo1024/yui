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
}
