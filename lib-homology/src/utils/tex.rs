//! LaTeX rendering for graded R-modules and graded objects.

use itertools::Itertools;
use yui_core::{Ring, RingOps, TeX};

use std::fmt::Display;

use crate::utils::format::make_rmod_str;
use crate::utils::{ToSeqString, ToTableString};

/// LaTeX form of [`ToSeqString`]: the same indices, entries rendered by [`TeX`].
pub trait ToTexSeq<I: Display>: ToSeqString<I> {
    fn tex_entry_at(&self, i: &I) -> String;

    fn tex_label(&self) -> String {
        self.label()
    }

    fn tex_seq(&self, caption: &str) -> String {
        yui_core::tex_table(caption, self.tex_label(), [""], self.indices(), |_, i| {
            self.tex_entry_at(i)
        }, true, true)
    }
}

/// LaTeX form of [`ToTableString`]: the same index ranges, so the two renderings agree,
/// with entries rendered by [`TeX`] instead of [`Display`].
pub trait ToTexTable<I: Display>: ToTableString<I> {
    fn tex_entry_at(&self, i: &I, j: &I) -> String;

    fn tex_labels(&self) -> (String, String) {
        self.labels()
    }

    fn tex_table(&self, caption: &str) -> String {
        let (label0, label1) = self.tex_labels();
        let (ind0, ind1) = self.indices();
        let head = format!("{label1} \\backslash {label0}");

        yui_core::tex_table(caption, head, ind1.into_iter().rev(), ind0, |j, i| {
            self.tex_entry_at(i, j)
        }, true, true)
    }
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
