use std::fmt::Display;
use std::collections::HashMap;
use std::ops::Index;
use itertools::Itertools;
use yui_core::{Ring, RingOps, EucRing, EucRingOps};
use yui_matrix::sparse::{SpMat, SpVec};
use yui_lin_comb::{FreeGen, LinComb};
use crate::utils::HomologyCalc;
use crate::{RModStr, GridItr, GridIdx, GenericRModGrid, Grid, ChainComplex, HomologyComputable, GenericHomology, GenericRModStr};
use crate::fmt::FmtList;

pub struct FreeRModStr<X, R>
where 
    X: FreeGen,
    R: Ring, for<'x> &'x R: RingOps<R>, 
{
    generators: Vec<X>,
    tors: Vec<R>, // always empty
}

impl<X, R> FreeRModStr<X, R>
where 
    X: FreeGen,
    R: Ring, for<'x> &'x R: RingOps<R>, 
{
    pub fn new(generators: Vec<X>) -> Self {
        Self { generators, tors: vec![] }
    }

    pub fn generators(&self) -> &Vec<X> {
        &self.generators
    }

    pub fn vectorize_x(&self, x: X) -> SpVec<R> {
        self.vectorize(&LinComb::from(x))
    }

    pub fn vectorize(&self, z: &LinComb<X, R>) -> SpVec<R> {
        let n = self.generators.len();

        SpVec::generate(n, |set| { 
            for (x, a) in z.iter() { 
                // MEMO non-effective
                let Some(i) = self.generators.iter().position(|z| x == z) else { 
                    continue
                };
                set(i, a.clone())
            }
        })
    }

    pub fn make_matrix<F>(&self, target: &Self, f: F) -> SpMat<R>
    where 
        F: Fn(&X) -> Vec<(X, R)> 
    {
        let (m, n) = (target.generators.len(), self.generators.len());

        let t_ind = target.generators.iter()
                .enumerate()
                .map(|(i, y)| (y, i))
                .collect::<HashMap<_, _>>();

        SpMat::generate((m, n), |set|
            for (j, x) in self.generators.iter().enumerate() {
                let ys = f(x);
                for (y, a) in ys {
                    if !t_ind.contains_key(&y) { continue }
                    let i = t_ind[&y];
                    set(i, j, a);
                }
            }
        )
    }
}

impl<X, R> RModStr for FreeRModStr<X, R>
where 
    X: FreeGen,
    R: Ring, for<'x> &'x R: RingOps<R>, 
{
    type R = R;

    fn zero() -> Self {
        Self::new(vec![])
    }

    fn rank(&self) -> usize {
        self.generators.len()
    }

    fn tors(&self) -> &Vec<Self::R> {
        &self.tors
    }
}

impl<X, R> Display for FreeRModStr<X, R>
where 
    X: FreeGen,
    R: Ring, for<'x> &'x R: RingOps<R>, 
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.fmt_default(f)
    }
}

pub struct FreeChainComplex<X, R, I>
where 
    X: FreeGen,
    R: Ring, for<'x> &'x R: RingOps<R>,
    I: GridItr, I::Item: GridIdx
{
    grid: GenericRModGrid<FreeRModStr<X, R>, I>,
    d_degree: I::Item,
    d_map: Box<dyn Fn(&X) -> Vec<(X, R)>>
}

impl<X, R, I> FreeChainComplex<X, R, I>
where 
    X: FreeGen,
    R: Ring, for<'x> &'x R: RingOps<R>,
    I: GridItr, I::Item: GridIdx
{
    pub fn new<F1, F2>(range: I, d_degree: I::Item, gens: F1, d_maps: F2) -> Self
    where F1: Fn(I::Item) -> Vec<X>, F2: Fn(&X) -> Vec<(X, R)> + 'static { 
        let grid = GenericRModGrid::new(range, |i| { 
            FreeRModStr::new(gens(i))
        });
        let d_maps = Box::new(d_maps);
        Self { grid, d_degree, d_map: d_maps }
    }

    pub fn differentiate_x(&self, x: &X) -> Vec<(X, R)> {
        (self.d_map)(x)
    }

    pub fn differetiate(&self, z: &LinComb<X, R>) -> LinComb<X, R> { 
        z.apply(|x| self.differentiate_x(x))
    }

    pub fn desc_d_gens_at(&self, i: I::Item) -> String {
        format!("C[{i}]: {}", &self[i]) + &self[i].generators.iter().map(|x| { 
            let dx = LinComb::from_iter(self.differentiate_x(x));
            format!("\n  {x} -> {dx}")
        }).join("")
    }

    pub fn desc_d_gens(&self) -> String { 
        self.indices().map(|i| self.desc_d_gens_at(i)).join("\n\n")
    }

    pub fn print_d_gens(&self) {
        println!("{}", self.desc_d_gens());
    }
}

impl<X, R, I> Display for FreeChainComplex<X, R, I>
where 
    X: FreeGen,
    R: Ring, for<'x> &'x R: RingOps<R>,
    I: GridItr, I::Item: GridIdx
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.fmt_list(f)
    }
}

impl<X, R, I> Index<I::Item> for FreeChainComplex<X, R, I>
where 
    X: FreeGen,
    R: Ring, for<'x> &'x R: RingOps<R>,
    I: GridItr, I::Item: GridIdx
{
    type Output = FreeRModStr<X, R>;

    fn index(&self, index: I::Item) -> &Self::Output {
        &self.grid[index]
    }
}

impl<X, R, I> Grid for FreeChainComplex<X, R, I>
where 
    X: FreeGen,
    R: Ring, for<'x> &'x R: RingOps<R>,
    I: GridItr, I::Item: GridIdx
{
    type Idx = I::Item;
    type IdxIter = I;
    type Output = FreeRModStr<X, R>;

    fn contains_idx(&self, k: Self::Idx) -> bool {
        self.grid.contains_idx(k)
    }

    fn indices(&self) -> Self::IdxIter {
        self.grid.indices()
    }

    fn get(&self, i: Self::Idx) -> Option<&Self::Output> {
        self.grid.get(i)
    }
}

impl<X, R, I> ChainComplex for FreeChainComplex<X, R, I>
where 
    X: FreeGen,
    R: Ring, for<'x> &'x R: RingOps<R>,
    I: GridItr, I::Item: GridIdx
{
    fn d_degree(&self) -> Self::Idx {
        self.d_degree
    }

    fn d_matrix(&self, k: Self::Idx) -> SpMat<Self::R> {
        let c0 = &self[k];
        let c1 = &self[k + self.d_degree];
        c0.make_matrix(c1, &self.d_map)
    }
}

impl<X, R, I> HomologyComputable for FreeChainComplex<X, R, I>
where 
    X: FreeGen,
    R: EucRing, for<'x> &'x R: EucRingOps<R>,
    I: GridItr,
    I::Item: GridIdx
{
    type Homology = GenericHomology<R, I>;
    type HomologySummand = GenericRModStr<R>;

    fn homology_at(&self, i: Self::Idx) -> Self::HomologySummand {
        let d1 = self.d_matrix(i - self.d_degree());
        let d2 = self.d_matrix(i);
        HomologyCalc::calculate(d1, d2)
    }

    fn homology(&self) -> Self::Homology {
        let range = self.indices();
        let grid = GenericRModGrid::new(range, |i| {
            self.homology_at(i)
        });
        Self::Homology::new(grid)
    }
}
