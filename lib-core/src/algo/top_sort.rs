use std::collections::VecDeque;
use ahash::AHashMap;
use num_traits::Zero;

pub fn top_sort<Itr>(itr: Itr) -> Result<Vec<usize>, String> 
where Itr: IntoIterator<Item = (usize, Vec<usize>)> {
    let data: AHashMap<usize, Vec<usize>> = itr.into_iter().collect();
    if data.is_empty() { 
        return Ok(vec![])
    }

    let mut weight: AHashMap<_, _> = data.keys().map(|&i| (i, 0)).collect();

    for s in data.values() { 
        for &i in s { 
            let Some(w) = weight.get_mut(&i) else { 
                Err(format!("Vertex {i} not in key."))?
            };
            *w += 1;
        }
    }

    let mut res = vec![];
    let mut queue = VecDeque::from_iter(
        weight.iter().filter_map(|(&i, w)| 
            if w.is_zero() { Some(i) } else { None }
        )
    );

    if queue.is_empty() {
        Err("Input is cyclic.")?
    }

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
    
    if res.len() < data.len() { 
        Err("Input contains cycle.")?
    }

    Ok(res)
}

#[cfg(test)]
mod tests {
    use super::top_sort;
 
    #[test]
    fn valid() { 
        // example from https://en.wikipedia.org/wiki/Topological_sorting
        let tree = vec![
            (5,  vec![11]),
            (7,  vec![11, 8]),
            (3,  vec![8, 10]),
            (11, vec![2, 9, 10]),
            (8,  vec![9]),
            (2,  vec![]),
            (9,  vec![]),
            (10, vec![])
        ];
        let res = top_sort(tree.clone());
        assert!(res.is_ok());

        let res = res.unwrap();
        assert!(tree.iter().all(|(i, list)| {
            list.iter().all(|j| {
                res.iter().position(|i1| i1 == i).unwrap()
                < res.iter().position(|j1| j1 == j).unwrap()
            })
        }));
    }

    #[test]
    fn missing_vertex() { 
        let tree = vec![
            (0, vec![1])
        ];
        let res = top_sort(tree);
        assert!(res.is_err());
    }

    #[test]
    fn cycle() { 
        let tree = vec![
            (0, vec![1]),
            (1, vec![0]),
        ];

        let res = top_sort(tree);
        assert!(res.is_err());
    }

    #[test]
    fn contains_cycle() { 
        let tree = vec![
            (0, vec![1]),
            (1, vec![2]),
            (2, vec![1]),
        ];

        let res = top_sort(tree);
        assert!(res.is_err());
    }
}