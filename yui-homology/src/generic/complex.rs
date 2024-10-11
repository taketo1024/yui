use yui::{Ring, RingOps};
use yui_matrix::sparse::SpMat;
use yui_matrix::MatTrait;
use super::gen::EnumGen;
use super::GenericSummand;

pub type GenericChainComplexBase<I, R> = ChainComplexBase<I, EnumGen<I>, R>;
pub type GenericChainComplex<R> = GenericChainComplexBase<isize, R>;

use crate::{Grid, GridDeg, GridTrait, ChainComplexBase};

impl<I, R> GenericChainComplexBase<I, R> 
where I: GridDeg, R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn generate<It, F>(support: It, d_deg: I, d_matrix_map: F) -> Self
    where 
        It: IntoIterator<Item = I>, 
        F: FnMut(I) -> SpMat<R>
    {
        let d_matrices = Grid::generate(support, d_matrix_map);

        let summands = Grid::generate(d_matrices.support(), |i| {
            let r = d_matrices[i].ncols();
            GenericSummand::generate_free(i, r)
        });

        Self::new(
            summands.clone(), d_deg, 
            move |i, z| {
                let d = &d_matrices[i];
                let v = summands[i].vectorize(z);
                let dv = d * v;
                let dz = summands[i + d_deg].as_chain(&dv);
                dz
            }
        )
    }
}

impl<R> GenericChainComplex<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    pub fn one() -> GenericChainComplex<R> {
        GenericChainComplex::generate(0..=0, -1, |i|
            match i { 
                0 => SpMat::from_dense_data((0, 1), []),
                _ => SpMat::zero((0, 0))
            }
        )
    }

    pub fn one_one(r: R) -> GenericChainComplex<R> {
        GenericChainComplex::generate(0..=1, -1, |i|
            match i { 
                0 => SpMat::from_dense_data((0, 1), []),
                1 => SpMat::from_dense_data((1, 1), [r.clone()]),
                _ => SpMat::zero((0, 0))
            }
        )
    }

    pub fn two_one(r1: R, r2: R) -> GenericChainComplex<R> {
        GenericChainComplex::generate(0..=1, -1, |i|
            match i { 
                0 => SpMat::from_dense_data((0, 1), []),
                1 => SpMat::from_dense_data((1, 2), [r1.clone(), r2.clone()]),
                _ => SpMat::zero((0, 0))
            }
        )
    }

    pub fn one_two(r1: R, r2: R) -> GenericChainComplex<R> {
        GenericChainComplex::generate(0..=1, -1, |i|
            match i { 
                0 => SpMat::from_dense_data((0, 2), []),
                1 => SpMat::from_dense_data((2, 1), [r1.clone(), r2.clone()]),
                _ => SpMat::zero((0, 0))
            }
        )
    }

    pub fn d3() -> GenericChainComplex<R> {
        GenericChainComplex::generate(0..=3, -1, |i|
            match i { 
                0 => Self::mat((0, 4), []),
                1 => Self::mat((4, 6), [-1, -1, 0, -1, 0, 0, 1, 0, -1, 0, -1, 0, 0, 1, 1, 0, 0, -1, 0, 0, 0, 1, 1, 1] ),
                2 => Self::mat((6, 4), [1, 1, 0, 0, -1, 0, 1, 0, 1, 0, 0, 1, 0, -1, -1, 0, 0, 1, 0, -1, 0, 0, 1, 1] ),
                3 => Self::mat((4, 1), [-1, 1, -1, 1]),
                _ => SpMat::zero((0, 0))
            }
        )
    }

    pub fn s2() -> GenericChainComplex<R>
     {
        GenericChainComplex::generate(0..=2, -1, |i|
            match i { 
                0 => Self::mat((0, 4), []),
                1 => Self::mat((4, 6), [-1, -1, 0, -1, 0, 0, 1, 0, -1, 0, -1, 0, 0, 1, 1, 0, 0, -1, 0, 0, 0, 1, 1, 1]),
                2 => Self::mat((6, 4), [1, 1, 0, 0, -1, 0, 1, 0, 1, 0, 0, 1, 0, -1, -1, 0, 0, 1, 0, -1, 0, 0, 1, 1] ),
                _ => SpMat::zero((0, 0))
            }
        )
    }

    pub fn t2() -> GenericChainComplex<R>
     {
        GenericChainComplex::generate(0..=2, -1, |i|
            match i { 
                0 => Self::mat((0, 9), []),
                1 => Self::mat((9, 27), [-1, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 1, 0, -1, 0, 0, 0, 0, 0, 0, -1, -1, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 1, -1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, -1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1]),
                2 => Self::mat((27, 18), [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1]),
                _ => SpMat::zero((0, 0))
            }
        )
    }

    pub fn rp2() -> GenericChainComplex<R>
     { 
        GenericChainComplex::generate(0..=2, -1, |i|
            match i { 
                0 => Self::mat((0, 6), []),
                1 => Self::mat((6, 15), [-1, -1, 0, 0, 0, 0, 0, -1, -1, 0, -1, 0, 0, 0, 0, 1, 0, -1, -1, 0, -1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 1, 1, 0, 1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 1, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1] ),
                2 => Self::mat((15, 10), [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1] ),
                _ => SpMat::zero((0, 0))
            }
        )
    }

    fn mat<I>(shape: (usize, usize), entries: I) -> SpMat<R>
    where I: IntoIterator<Item = i32> { 
        SpMat::from_dense_data(shape, entries.into_iter().map(|a| 
            match a { 
                0 => R::zero(),
                1 => R::one(),
                -1 => -R::one(),
                _ => panic!()
            }
        ))
    }
}

#[cfg(test)]
pub(crate) mod tests { 
    use num_traits::Zero;

    use crate::{ChainComplexTrait, SummandTrait};

    use super::*;

    #[test]
    fn zero() { 
        let c = GenericChainComplex::<i32>::zero();
        assert_eq!(c[0].rank(), 0);
        
        c.check_d_all();

        let h = c.homology();
        assert!(h[0].is_zero());
    }

    #[test]
    fn single() { 
        let c = GenericChainComplex::<i32>::one();
        assert_eq!(c[0].rank(), 1);

        c.check_d_all();

        let h = c.homology();
        
        assert_eq!(h[0].rank(), 1);
        assert!( h[0].is_free());
        assert!(!h[0].is_zero());
    }

    #[test]
    fn one_to_one() { 
        let c = GenericChainComplex::<i32>::one_one(1);
        let h = c.homology();

        assert!(h[0].is_zero());
        assert!(h[1].is_zero());
    }

    #[test]
    fn two_to_one() { 
        let c = GenericChainComplex::<i32>::two_one(1, -1);
        let h = c.homology();

        assert!(h[0].is_zero());
        assert_eq!(h[1].rank(), 1);
        assert!(h[1].is_free());
    }

    #[test]
    fn one_to_two() { 
        let c = GenericChainComplex::<i32>::one_two(1, -1);
        let h = c.homology();

        assert_eq!(h[0].rank(), 1);
        assert!(h[0].is_free());
        assert!(h[1].is_zero());
    }

    #[test]
    fn torsion() { 
        let c = GenericChainComplex::<i32>::one_one(2);
        let h = c.homology();

        assert_eq!(h[0].rank(), 0);
        assert_eq!(h[0].tors(), &vec![2]);
        assert!(!h[0].is_free());
    }

    #[test]
    fn d3() {
        let c = GenericChainComplex::<i32>::d3();

        assert_eq!(c[0].rank(), 4);
        assert_eq!(c[1].rank(), 6);
        assert_eq!(c[2].rank(), 4);
        assert_eq!(c[3].rank(), 1);

        c.check_d_all();

        let h = c.homology();

        assert_eq!(h[0].rank(), 1);
        assert!(h[0].is_free());

        assert_eq!(h[1].rank(), 0);
        assert!(h[1].is_free());

        assert_eq!(h[2].rank(), 0);
        assert!(h[2].is_free());

        assert_eq!(h[3].rank(), 0);
        assert!(h[3].is_free());
    }

    #[test]
    fn s2() {
        let c = GenericChainComplex::<i32>::s2();

        assert_eq!(c[0].rank(), 4);
        assert_eq!(c[1].rank(), 6);
        assert_eq!(c[2].rank(), 4);
        assert_eq!(c[3].rank(), 0);

        c.check_d_all();

        let h = c.homology();

        assert_eq!(h[0].rank(), 1);
        assert!(h[0].is_free());

        assert_eq!(h[1].rank(), 0);
        assert!(h[1].is_free());

        assert_eq!(h[2].rank(), 1);
        assert!(h[2].is_free());
    }

    #[test]
    fn t2() {
        let c = GenericChainComplex::<i32>::t2();

        assert_eq!(c[0].rank(), 9);
        assert_eq!(c[1].rank(), 27);
        assert_eq!(c[2].rank(), 18);
        assert_eq!(c[3].rank(), 0);

        c.check_d_all();

        let h = c.homology();

        assert_eq!(h[0].rank(), 1);
        assert!(h[0].is_free());

        assert_eq!(h[1].rank(), 2);
        assert!(h[1].is_free());

        assert_eq!(h[2].rank(), 1);
        assert!(h[2].is_free());
    }

    #[test]
    fn rp2() {
        let c = GenericChainComplex::<i32>::rp2();

        assert_eq!(c[0].rank(), 6);
        assert_eq!(c[1].rank(), 15);
        assert_eq!(c[2].rank(), 10);
        assert_eq!(c[3].rank(), 0);

        c.check_d_all();
        
        let h = c.homology();

        assert_eq!(h[0].rank(), 1);
        assert!(h[0].is_free());

        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[1].tors(), &vec![2]);
        assert!(!h[1].is_free());

        assert_eq!(h[2].rank(), 0);
        assert!(h[2].is_free());
    }

    #[test]
    fn s2_gens() {
        let c = GenericChainComplex::<i32>::s2();
        let h = c.homology();

        let z = h[2].gen_chain(0);
        let dz = c.d(2, &z);

        assert!(!z.is_zero());
        assert!(dz.is_zero());
    }

    #[test]
    fn t2_gens() {
        let c = GenericChainComplex::<i32>::t2();
        let h = c.homology();

        let z = h[2].gen_chain(0);
        let dz = c.d(2, &z);

        assert!(!z.is_zero());
        assert!(dz.is_zero());

        let a = h[1].gen_chain(0);
        let b = h[1].gen_chain(1);
        let da = c.d(1, &a);
        let db = c.d(1, &b);

        assert!(!a.is_zero());
        assert!(!b.is_zero());
        assert!(da.is_zero());
        assert!(db.is_zero());
    }

    #[test]
    fn rp2_gens() {
        let c = GenericChainComplex::<i32>::rp2();
        let h = c.homology();

        let z = h[1].gen_chain(0);
        let dz = c.d(1, &z);

        assert!(!z.is_zero());
        assert!(dz.is_zero());

        assert!(!h[1].vectorize_euc(&z).is_zero());
        assert!(h[1].vectorize_euc(&(z * 2)).is_zero()); // order 2
    }

    #[test]
    fn d3_red() {
        let c = GenericChainComplex::<i32>::d3().reduced();

        assert_eq!(c[0].rank(), 1);
        assert_eq!(c[1].rank(), 0);
        assert_eq!(c[2].rank(), 0);
        assert_eq!(c[3].rank(), 0);

        c.check_d_all();

        let h = c.homology();

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[2].rank(), 0);
        assert_eq!(h[3].rank(), 0);

        let z = h[0].gen_chain(0);
        assert!(!z.is_zero());
        assert!(c.d(0, &z).is_zero());
    }

    #[test]
    fn s2_red() {
        let c = GenericChainComplex::<i32>::s2().reduced();

        assert_eq!(c[0].rank(), 1);
        assert_eq!(c[1].rank(), 0);
        assert_eq!(c[2].rank(), 1);

        c.check_d_all();

        let h = c.homology();

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[2].rank(), 1);

        for i in 0..=2 { 
            for j in 0..h[i].rank() { 
                let z = h[i].gen_chain(j);
                assert!(!z.is_zero());
                assert!(c.d(i, &z).is_zero());
            }
        }
    }

    #[test]
    fn t2_red() {
        let c = GenericChainComplex::<i32>::t2().reduced();

        assert_eq!(c[0].rank(), 1);
        assert_eq!(c[1].rank(), 2);
        assert_eq!(c[2].rank(), 1);

        c.check_d_all();

        let h = c.homology();

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[1].rank(), 2);
        assert_eq!(h[2].rank(), 1);

        for i in 0..=2 { 
            for j in 0..h[i].rank() { 
                let z = h[i].gen_chain(j);
                assert!(!z.is_zero());
                assert!(c.d(i, &z).is_zero());
            }
        }
    }

    #[test]
    fn rp2_red() {
        let c = GenericChainComplex::<i32>::rp2().reduced();

        assert_eq!(c[0].rank(), 1);
        assert_eq!(c[1].rank(), 1);
        assert_eq!(c[2].rank(), 1);

        c.check_d_all();
        
        let h = c.homology();

        assert_eq!(h[0].rank(), 1);
        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[1].tors(), &vec![2]);
        assert_eq!(h[2].rank(), 0);

        for i in 0..=2 { 
            for j in 0..h[i].rank() { 
                let z = h[i].gen_chain(j);
                assert!(!z.is_zero());
                assert!(c.d(i, &z).is_zero());
            }
        }
    }
}