use std::ops::RangeInclusive;

pub fn range_of<Idx, Itr>(itr: Itr) -> RangeInclusive<Idx>
where Idx: Ord + Default + Copy, Itr: IntoIterator<Item = Idx> { 
    let (min, max) = itr.into_iter().fold(None, |res, i| { 
        if let Some((min, max)) = res { 
            if i < min { 
                Some((i, max))
            } else if max < i { 
                Some((min, i))
            } else { 
                Some((min, max))
            }
        } else {
            Some((i, i))
        }
    }).unwrap_or((Idx::default(), Idx::default()));

    min ..= max
}