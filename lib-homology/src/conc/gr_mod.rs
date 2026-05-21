use std::ops::{Index, RangeInclusive};

use ahash::AHashMap;
use itertools::Itertools;
use yui_core::lc::LcKey;
use yui_core::{Ring, RingOps};

use crate::utils::{ToSeqString, ToTableString};
use crate::{AddInd, Summand, isize2, usize2, isize3, usize3};

pub type GrMod1<X, R> = GrMod<isize,  X, R>;
pub type GrMod2<X, R> = GrMod<isize2, X, R>;
pub type GrMod3<X, R> = GrMod<isize3, X, R>;

/// An `I`-graded R-module: a sparse map from indices `I` to free / finitely
/// generated R-modules `Summand<X, R>`. Missing entries are treated as the
/// zero module via the `default` field.
#[derive(Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct GrMod<I, X, R>
where
    I: AddInd,
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>,
{
    data: AHashMap<I, Summand<X, R>>,
    #[cfg_attr(feature = "serde", serde(skip))]
    zero_summand: Summand<X, R>,
}

impl<I, X, R> GrMod<I, X, R>
where
    I: AddInd,
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>,
{
    fn new(data: impl IntoIterator<Item = (I, Summand<X, R>)>) -> Self {
        Self { data: data.into_iter().collect(), zero_summand: Summand::zero() }
    }

    pub fn generate<It, F>(support: It, mut e_map: F) -> Self
    where
        It: IntoIterator<Item = I>,
        F: FnMut(I) -> Summand<X, R>,
    {
        Self::new(support.into_iter().map(|i| (i, e_map(i))))
    }

    pub fn generate_filtered<It, F>(support: It, mut e_map: F) -> Self
    where
        It: IntoIterator<Item = I>,
        F: FnMut(I) -> Option<Summand<X, R>>,
    {
        Self::new(support.into_iter().filter_map(|i| e_map(i).map(|e| (i, e))))
    }

    pub fn support(&self) -> impl Iterator<Item = &I> + '_ {
        self.data.keys()
    }

    pub fn is_supported(&self, i: I) -> bool {
        self.data.contains_key(&i)
    }

    pub fn get(&self, i: I) -> Option<&Summand<X, R>> {
        self.data.get(&i)
    }

    pub fn zero_summand(&self) -> &Summand<X, R> {
        &self.zero_summand
    }

    pub fn iter(&self) -> impl Iterator<Item = (I, &Summand<X, R>)> {
        self.data.iter().map(|(&i, e)| (i, e))
    }

    pub fn total_rank(&self) -> usize {
        self.support().map(|&i| self[i].rank()).sum()
    }
}

impl<X, R> GrMod1<X, R>
where X: LcKey, R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn truncated(&self, range: RangeInclusive<isize>) -> Self {
        Self::new(
            self.data.iter().filter_map(|(&i, e)| range.contains(&i).then_some((i, e.clone())))
        )
    }
}

impl<I, X, R> Default for GrMod<I, X, R>
where
    I: AddInd,
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>,
{
    fn default() -> Self {
        Self::new(AHashMap::default())
    }
}

impl<I, X, R> Index<I> for GrMod<I, X, R>
where
    I: AddInd,
    X: LcKey,
    R: Ring, for<'x> &'x R: RingOps<R>,
{
    type Output = Summand<X, R>;
    fn index(&self, i: I) -> &Self::Output {
        self.data.get(&i).unwrap_or(&self.zero_summand)
    }
}

macro_rules! impl_index {
    ($t:ident, $t2:ident, $t3:ident) => {
        impl<X, R> Index<($t, $t)> for GrMod<$t2, X, R>
        where X: LcKey, R: Ring, for<'x> &'x R: RingOps<R> {
            type Output = Summand<X, R>;
            fn index(&self, i: ($t, $t)) -> &Self::Output {
                &self[$t2::from(i)]
            }
        }

        impl<X, R> Index<($t, $t, $t)> for GrMod<$t3, X, R>
        where X: LcKey, R: Ring, for<'x> &'x R: RingOps<R> {
            type Output = Summand<X, R>;
            fn index(&self, i: ($t, $t, $t)) -> &Self::Output {
                &self[$t3::from(i)]
            }
        }
    };
}

impl_index!(isize, isize2, isize3);
impl_index!(usize, usize2, usize3);

impl<X, R> ToSeqString<isize> for GrMod1<X, R>
where X: LcKey, R: Ring, for<'x> &'x R: RingOps<R> {
    fn label(&self) -> String {
        "i".to_string()
    }

    fn indices(&self) -> Vec<isize> {
        self.support().sorted().cloned().collect()
    }

    fn entry_at(&self, i: &isize) -> String {
        self.get(*i).map(|e| e.to_string()).unwrap_or_else(|| ".".to_string())
    }
}

impl<X, R> ToTableString<isize> for GrMod2<X, R>
where X: LcKey, R: Ring, for<'x> &'x R: RingOps<R> {
    fn labels(&self) -> (String, String) {
        ("i".to_string(), "j".to_string())
    }

    fn indices(&self) -> (Vec<isize>, Vec<isize>) {
        let is = self.support().map(|&isize2(i, _)| i).unique().sorted().collect();
        let js = self.support().map(|&isize2(_, j)| j).unique().sorted().collect();
        (is, js)
    }

    fn entry_at(&self, i: &isize, j: &isize) -> String {
        self.get(isize2(*i, *j)).map(|e| e.to_string()).unwrap_or_else(|| ".".to_string())
    }
}

#[cfg(feature = "tex")]
mod tex_impl {
    use super::*;
    use yui_core::{TeX, tex_table};
    use crate::utils::tex::TeXTable;

    macro_rules! impl_tex_table {
        ($t:ident) => {
            impl<X, R> TeXTable<$t> for GrMod<$t, X, R>
            where X: LcKey, R: Ring + TeX, for<'x> &'x R: RingOps<R> {
                fn tex_table(&self, caption: &str, head: &str) -> String {
                    let cols = self.support().map(|&$t(i, _)| i).unique().sorted();
                    let rows = self.support().map(|&$t(_, j)| j).unique().sorted().rev();

                    tex_table(caption, head, rows, cols, |&j, &i| {
                        self.get($t(i, j))
                            .map(|e| e.tex_string())
                            .unwrap_or_else(|| ".".to_string())
                    }, true, false)
                }
            }
        };
    }

    impl_tex_table!(isize2);
}
