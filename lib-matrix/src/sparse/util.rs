pub(crate) fn perm_for_indices<'a, I>(n: usize, indices: I) -> sprs::PermOwned
where I: IntoIterator<Item = &'a usize> { 
    use std::collections::BTreeSet;

    let mut set: BTreeSet<_> = (0..n).collect();
    let mut vec: Vec<usize> = vec![];

    for &i in indices { 
        assert!(i < n);
        vec.push(i);
        set.remove(&i);
    }
    for i in set { 
        vec.push(i);
    }

    let mut inv = vec![0; n];
    for (i, j) in vec.into_iter().enumerate() {
        inv[j] = i;
    }

    sprs::PermOwned::new(inv)
}