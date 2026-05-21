use itertools::Itertools;
use yui_core::{Ring, RingOps, TeX};

use crate::utils::format::make_rmod_str;

pub trait TeXTable<I> {
    fn tex_table(&self, caption: &str, head: &str) -> String;
}

pub fn tex_rmod_str<R>(rank: usize, tors: &[R]) -> String
where R: Ring + TeX, for<'x> &'x R: RingOps<R> {
    make_rmod_str(
        R::tex_math_symbol(),
        rank,
        &tors.iter().map(|t| t.tex_string()).collect_vec(),
        |a| format!("^{{{a}}}"),
        "\\oplus"
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tex() {
        let s = tex_rmod_str(2, &[2,2,3]);
        assert_eq!(s, "\\mathbb{Z}^{2} \\oplus (\\mathbb{Z}/2)^{2} \\oplus (\\mathbb{Z}/3)");
    }
}
