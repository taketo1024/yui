pub trait MatTrait {
    fn shape(&self) -> (usize, usize);
    fn nrows(&self) -> usize { self.shape().0 }
    fn ncols(&self) -> usize { self.shape().1 }
    fn is_square(&self) -> bool { 
        let (m, n) = self.shape();
        m == n
    }
}