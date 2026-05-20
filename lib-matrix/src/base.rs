pub trait MatTrait {
    fn shape(&self) -> (usize, usize);
    fn n_rows(&self) -> usize { self.shape().0 }
    fn n_cols(&self) -> usize { self.shape().1 }
    fn is_square(&self) -> bool {
        let (m, n) = self.shape();
        m == n
    }
}