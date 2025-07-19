//! Simple yet efficient implementation of a k-d tree.
//!
//! A k-d tree is a binary search tree where data in each node is a k-dimensional point in space.
//! Once a tree has been constructed by the [`create_kd_tree()`](create_kd_tree())) function,
//! it allows to quickly find a nearest neighbor of a k-dimensional query
//! (performed by [`find_nearest()`](find_nearest())) as well as to locate all neighbors
//! within a given distance from a query. For the latter use (performed by [`find_within()`](find_within())).
//!
//! ```rust
//! use bioshell_datastructures::euclidean_distance_squared;
//! use bioshell_datastructures::kd_tree::{create_kd_tree, find_nearest, find_within};
//! let mut points = [[0.1, 0.2], [0.2, 0.2], [1.1, 1.2], [2.2, 2.2]];
//! let root = create_kd_tree(&mut points, 2).unwrap();
//! let query = [0.3, 0.3];
//! let (d, e) = find_nearest(&root, &query, 2, euclidean_distance_squared);
//! assert!((d - 0.02).abs() < 0.000001);
//! let neighbors = find_within(&root, &query, 2, 0.1, euclidean_distance_squared);
//! assert_eq!(neighbors.len(),2);
//! ```

use std::ops::Index;

use crate::BinaryTreeNode;


pub struct KdTreeData<T> {
    /// k-dimensional data element for the branching point
    pub value: T,
    /// level of this node in the tree; level of the root is 0
    pub level: usize,
    /// which coordinate was used to split the data at this node
    pub split_coordinate: usize
}


pub fn create_kd_tree<T>(data: &mut [T], dimensionality: usize) -> Option<Box<BinaryTreeNode<KdTreeData<T>>>>
    where T: Index<usize, Output = f64> + std::fmt::Debug, T:Clone {

    fn create_kd_tree_rec<T>(data: &mut [T], tree_depth: usize, dimensionality: usize, next_id: u32) -> Option<Box<BinaryTreeNode<KdTreeData<T>>>>
        where T: Index<usize, Output = f64> + std::fmt::Debug, T:Clone {

        if data.len()==0 { return None; }
        if data.len() == 1 {
            let kdpoint = KdTreeData{value: data[0].clone(), level: tree_depth+1, split_coordinate: 0 };
            let mut node = BinaryTreeNode::new(kdpoint);
            node.id = next_id;
            return Some(Box::new(node));
        }

        sort_along_dimension(data, tree_depth % dimensionality);
        let median = data.len() / 2;
        let root_data = KdTreeData{ value: data[median].clone(), level: tree_depth+1, split_coordinate: tree_depth % dimensionality };
        let mut root = BinaryTreeNode::new(root_data);
        root.id = next_id;
        let (left_data, mut right_data) = data.split_at_mut(median);
        right_data = &mut right_data[1..];
        root.left = create_kd_tree_rec(left_data, tree_depth + 1, dimensionality, next_id * 2);
        root.right = create_kd_tree_rec(right_data, tree_depth + 1, dimensionality, next_id * 2 + 1);

        return Some(Box::new(root));
    }

    return create_kd_tree_rec(data, 0, dimensionality, 1);
}


/// Finds the point nearest to a given query.
///
/// # Arguments
/// * `tree_root` - k-dimensional points is borrowed mutably since the points will be re-ordered during this call
/// * `query` - a query point
/// * `dimensionality` - the number of dimensions for each ot point
/// * `distance` - distance function
///
/// ```rust
/// use bioshell_datastructures::euclidean_distance_squared;
/// use bioshell_datastructures::kd_tree::{create_kd_tree, find_nearest};
/// let mut points = [[0.1, 0.2], [0.2, 0.2], [1.1, 1.2], [2.2, 2.2]];
/// let root = create_kd_tree(&mut points, 2).unwrap();
/// let query = [0.3, 0.3];
/// let (d, e) = find_nearest(&root, &query, 2, euclidean_distance_squared);
/// assert!((d - 0.02).abs() < 0.000001);
/// ```
pub fn find_nearest<'a, T, F>(tree_root: &'a Box<BinaryTreeNode<KdTreeData<T>>>, query: &T, dimensionality: usize, distance:F) -> (f64, &'a T)
    where T: Index<usize, Output = f64>, F: Fn(&T, &T, usize) -> f64 {

    let mut min_dist: f64 = distance(&tree_root.value.value, query, dimensionality);
    let mut min_elem: &T = &tree_root.value.value;
    let mut stack: Vec<&Box<BinaryTreeNode<KdTreeData<T>>>> = vec![];
    stack.push(tree_root);
    while !stack.is_empty() {
        let n = stack.pop().unwrap();
        let ne = &n.value.value;
        let d = distance(ne, query, dimensionality);
        if d < min_dist { (min_dist, min_elem) = (d, ne); }

        let k = n.value.split_coordinate;
        let d = (ne[k] - query[k]) * (ne[k] - query[k]);

        if d > min_dist {
            if ne[k] > query[k] {
                if let Some(left) = &n.left { stack.push(&left); }
            } else {
                if let Some(right) = &n.right { stack.push(&right); }
            }
        } else {
            if let Some(left) = &n.left { stack.push(&left); }
            if let Some(right) = &n.right { stack.push(&right); }
        }
    }

    return (min_dist, min_elem);
}

/// Finds all the point that are within a given radius from a given query
///
/// # Arguments
/// * `tree_root` - k-dimensional points is borrowed mutably since the points will be re-ordered during this call
/// * `query` - a query point
/// * `dimensionality` - the number of dimensions for each ot point
/// * `radius` - the search distance cutoff
/// * `distance` - distance function
///
/// ```rust
/// use bioshell_datastructures::kd_tree::{create_kd_tree, find_within};
/// use bioshell_datastructures::euclidean_distance_squared;
/// let mut points = [[0.1, 0.2], [0.2, 0.2], [1.1, 1.2], [2.2, 2.2]];
/// let root = create_kd_tree(&mut points, 2).unwrap();
/// let query = [0.3, 0.3];
/// let neighbors = find_within(&root, &query, 2, 0.1, euclidean_distance_squared);
/// assert_eq!(neighbors.len(), 2);
/// ```
pub fn find_within<'a, T, F>(tree_root: &'a Box<BinaryTreeNode<KdTreeData<T>>>, query: &T, dimensionality: usize, radius: f64, distance:F) -> Vec<&'a T>
    where T: Index<usize, Output = f64>, F: Fn(&T, &T, usize) -> f64 {

    let mut ret: Vec<&T> = vec![];

    let mut stack: Vec<&Box<BinaryTreeNode<KdTreeData<T>>>> = vec![];
    stack.push(tree_root);
    while !stack.is_empty() {
        let n = stack.pop().unwrap();
        let ne = &n.value.value;
        let d = distance(ne, query, dimensionality);
        if d < radius { ret.push(&n.value.value); }

        let k = n.value.split_coordinate;
        let d = (ne[k]-query[k])*(ne[k]-query[k]);
        if d > radius {
            if ne[k] > query[k] {
                if let Some(left) = &n.left { stack.push(&left); }
            } else {
                if let Some(right) = &n.right { stack.push(&right); }
            }
        } else {
            if let Some(left) = &n.left { stack.push(&left); }
            if let Some(right) = &n.right { stack.push(&right); }
        }
    }

    return ret;
}

fn sort_along_dimension<T>(data: &mut [T], idx:usize)  where T: Index<usize, Output = f64> {
    data.sort_by(|a,b|a[idx].partial_cmp(&b[idx]).unwrap());
}
