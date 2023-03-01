use std::collections::{HashMap, VecDeque};

use num_traits::Zero;


pub fn top_sort(data: HashMap<usize, Vec<usize>>) -> Vec<usize> {
    if data.is_empty() { 
        return vec![]
    }

    let mut weight: HashMap<_, _> = data.keys().map(|&i| (i, 0)).collect();

    for s in data.values() { 
        for &i in s { 
            let w = weight.get_mut(&i).unwrap();
            *w += 1;
        }
    }

    let mut res = vec![];
    let mut queue = VecDeque::from_iter(
        weight.iter().filter_map(|(&i, w)| 
            if w.is_zero() { Some(i) } else { None }
        )
    );

    assert!(!queue.is_empty());

    while let Some(i) = queue.pop_front() { 
        res.push(i);

        for &j in data[&i].iter() {
            let w = weight.get_mut(&j).unwrap();
            *w -= 1;

            if w.is_zero() { 
                queue.push_back(j);
            }
        }
    }
    
    assert_eq!(res.len(), data.len());

    res
}