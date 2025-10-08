use ahash::AHashSet;

use crate::{Ring, RingOps};

pub fn naive_det<R>(n: usize, matrix: &[R]) -> R 
where R: Ring, for<'x> &'x R: RingOps<R> {
    assert_eq!(matrix.len(), n * n);

    fn det_rec<R>(n: usize, matrix: &[R], used: AHashSet<usize>) -> R
    where R: Ring, for<'x> &'x R: RingOps<R> {
        let i = used.len();
        if i == n { 
            return R::one();
        }

        let mut e = -R::one();

        R::sum((0..n).flat_map(move |j| { 
            if used.contains(&j) { 
                return None
            }

            let a = &matrix[n * i + j];
            e = -&e;

            if a.is_zero() { 
                return None
            }

            let mut next = used.clone();
            next.insert(j);

            let d = &e * a * det_rec(n, matrix, next);
            
            Some(d)
        }))
    }

    det_rec(n, matrix, [].into())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_naive_det_1x1() {
        let matrix = [5];
        assert_eq!(naive_det(1, &matrix), 5);
    }

    #[test]
    fn test_naive_det_2x2() {
        let matrix = [
            1, 2, 
            3, 4
        ];
        assert_eq!(naive_det(2, &matrix), -2);
    }

    #[test]
    fn test_naive_det_3x3() {
        let matrix = [
            6, 1, 1,
            4, -2, 5,
            2, 8, 7
        ];
        assert_eq!(naive_det(3, &matrix), -306);
    }

    #[test]
    fn test_naive_det_4x4() {
        let matrix = [
            3, 2, 0, 1,
            4, 0, 1, 2,
            3, 0, 2, 1,
            9, 2, 3, 1
        ];
        assert_eq!(naive_det(4, &matrix), 24);
    }

    #[test]
    fn test_naive_det_5x5() {
        let matrix = [
            2, 0, 1, 3, 4,
            1, 2, 0, 1, 5,
            3, 1, 2, 1, 0,
            0, 2, 3, 2, 1,
            4, 1, 0, 2, 3
        ];
        assert_eq!(naive_det(5, &matrix), -150);
    }

    #[test]
    fn test_naive_det_identity() {
        // 3x3 identity matrix
        let matrix = [
            1, 0, 0,
            0, 1, 0,
            0, 0, 1
        ];
        assert_eq!(naive_det(3, &matrix), 1);
    }

    #[test]
    fn test_naive_det_zero_matrix() {
        let matrix = [0, 0, 0, 0];
        assert_eq!(naive_det(2, &matrix), 0);
    }

    #[test]
    fn test_naive_det_size_zero_matrix() {
        let matrix: [i32; 0] = [];
        assert_eq!(naive_det(0, &matrix), 1);
    }
}
