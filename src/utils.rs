pub fn ensure_len<T: Default>(v: &mut Vec<T>, n: usize) {
    if v.len() < n {
        v.resize_with(n, Default::default)
    }
}
