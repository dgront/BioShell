use std::ops::Index;
use num_traits::Float;

#[derive(Clone, Debug)]
pub struct KdTreeNode<T> {
    value: T,
    left: Option<Box<KdTreeNode<T>>>,
    right: Option<Box<KdTreeNode<T>>>,
    level: usize,
    split_coordinate: usize
}

pub fn create_kd_tree<T>(data: &mut [T], tree_depth: usize, dimensionality: usize) -> Option<Box<KdTreeNode<T>>>
    where <T as Index<usize>>::Output: Float, T: std::ops::Index<usize>, T:Clone {

    if data.len()==0 { return None; }
    if data.len() == 1 {
        let n = KdTreeNode{ value: data[0].clone(), left: None, right: None, level: tree_depth+1, split_coordinate: 0 };
        return Some(Box::new(n));
    }

    sort_along_dimension(data, tree_depth % dimensionality);
    let median = data.len() / 2;
    let mut root = KdTreeNode{ value: data[median].clone(), left: None, right: None, level: tree_depth+1, split_coordinate: tree_depth % dimensionality };
    let (left_data, right_data) = data.split_at_mut(median);
    root.left = create_kd_tree(left_data, tree_depth+1, dimensionality);
    root.right = create_kd_tree(right_data, tree_depth+1, dimensionality);

    return Some(Box::new(root));
}

fn sort_along_dimension<T>(data: &mut [T], idx:usize)
    where <T as Index<usize>>::Output: Float, T: std::ops::Index<usize> {
    data.sort_by(|a,b|b[idx].partial_cmp(&a[idx]).unwrap());
}