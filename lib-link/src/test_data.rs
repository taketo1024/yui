//! Hard-coded PD codes and braid words for knots up to 6 crossings,
//! used only in `cargo test` (and in downstream crates that opt into
//! the `test-utils` feature). Data comes from the KnotAtlas Take Home
//! Database / KnotInfo — its only use is to drive tests, so the
//! choice of chirality convention is self-contained.

use crate::{Braid, Link};

impl Link {
    pub fn test_data(name: &str) -> Result<Link, String> {
        let pd: &[[usize; 4]] = match name {
            "3_1" => &[[1,4,2,5],[3,6,4,1],[5,2,6,3]],
            "4_1" => &[[4,2,5,1],[8,6,1,5],[6,3,7,4],[2,7,3,8]],
            "5_1" => &[[1,6,2,7],[3,8,4,9],[5,10,6,1],[7,2,8,3],[9,4,10,5]],
            "5_2" => &[[1,4,2,5],[3,8,4,9],[5,10,6,1],[9,6,10,7],[7,2,8,3]],
            "6_1" => &[[1,4,2,5],[7,10,8,11],[3,9,4,8],[9,3,10,2],[5,12,6,1],[11,6,12,7]],
            "6_2" => &[[1,4,2,5],[5,10,6,11],[3,9,4,8],[9,3,10,2],[7,12,8,1],[11,6,12,7]],
            "6_3" => &[[4,2,5,1],[8,4,9,3],[12,9,1,10],[10,5,11,6],[6,11,7,12],[2,8,3,7]],
            "L2a1" => &[[4,1,3,2],[2,3,1,4]],
            _ => return Err(format!("no test data for `{name}`")),
        };
        Ok(Link::from_pd_code(pd.iter().copied()))
    }
}

impl Braid {
    pub fn test_data(name: &str) -> Result<Braid, String> {
        let word: &[i32] = match name {
            "3_1" => &[1, 1, 1],
            "4_1" => &[1, -2, 1, -2],
            "5_1" => &[1, 1, 1, 1, 1],
            "5_2" => &[1, 1, 1, 2, -1, 2],
            "6_1" => &[1, 1, 2, -1, -3, 2, -3],
            "6_2" => &[1, 1, 1, -2, 1, -2],
            "6_3" => &[1, 1, -2, 1, -2, -2],
            _ => return Err(format!("no test data for `{name}`")),
        };
        Ok(Braid::from_iter(word.iter().copied()))
    }
}
