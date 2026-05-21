//! LaTeX rendering for graded R-modules and graded objects, feature-gated
//! behind `tex`.

use itertools::Itertools;
use yui_core::{Ring, RingOps, TeX};

use crate::utils::format::make_rmod_str;

/// Render a bi-graded object as a LaTeX `tabular` table.
pub trait TeXTable<I> {
    fn tex_table(&self, caption: &str, head: &str) -> String;
}

/// LaTeX form of [`rmod_str`](crate::rmod_str). Renders to a math-mode string
/// like `\mathbb{Z}^{2} \oplus (\mathbb{Z}/2)^{2} \oplus (\mathbb{Z}/3)`.
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
