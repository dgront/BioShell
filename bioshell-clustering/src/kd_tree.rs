//! Simple yet efficient implementation of a k-d tree.
//!
//! A k-d tree is a binary search tree where data in each node is a k-dimensional point in space.
//!
//!
//! ```rust
//! use bioshell_clustering::{euclidean_distance_squared};
//! use bioshell_clustering::kd_tree::{create_kd_tree, find_nearest, find_within};
//! let mut points = [[0.1, 0.2], [0.2, 0.2], [1.1, 1.2], [2.2, 2.2]];
//! let root = create_kd_tree(&mut points, 2).unwrap();
//! let query = [0.3, 0.3];
//! let (d, e) = find_nearest(&root, &query, 2, euclidean_distance_squared);
//! assert!((d - 0.02).abs() < 0.000001);
//! let neighbors = find_within(&root, &query, 2, 0.1, euclidean_distance_squared);
//! assert_eq!(neighbors.len(),2);
//! ```

use std::ops::Index;

/// Represents a node of a K-dimensional tree.
#[derive(Clone)]
pub struct KdTreeNode<T> {
    /// index of the node
    id: u32,
    /// k-dimensional data element for the branching point
    value: T,
    /// left sub-tree holds all points located to the left of the branching point
    left: Option<Box<KdTreeNode<T>>>,
    /// right sub-tree holds all points located to the right of the branching point
    right: Option<Box<KdTreeNode<T>>>,
    /// level of this node in the tree; level of the root is 0
    level: usize,
    /// which coordinate was used to split the data at this node
    split_coordinate: usize
}

impl<T> KdTreeNode<T> {
    /// k-dimensional data element stored in this node
    ///
    /// This point has been used as the branching point at this node
    pub fn element(&self) -> &T { &self.value }
}

/// Creates a new k-d tree for a given array of k-dimensional points.
///
/// # Arguments
/// * `data` - k-dimensional points is borrowed mutably since the points will be re-ordered during this call
/// * `dimensionality` - the number of dimensions for each ot point of the generic type `T`
///
/// ```rust
/// use bioshell_clustering::kd_tree::create_kd_tree;
/// let mut points = [[0.1, 0.2], [0.2, 0.2], [1.1, 1.2], [2.2, 2.2]];
/// let root = create_kd_tree(&mut points, 2);
/// ```
///
pub fn create_kd_tree<T>(data: &mut [T], dimensionality: usize) -> Option<Box<KdTreeNode<T>>>
    where T: Index<usize, Output = f64> + std::fmt::Debug, T:Clone {

    fn create_kd_tree_rec<T>(data: &mut [T], tree_depth: usize, dimensionality: usize, next_id: u32) -> Option<Box<KdTreeNode<T>>>
        where T: Index<usize, Output = f64> + std::fmt::Debug, T:Clone {

        if data.len()==0 { return None; }
        if data.len() == 1 {
            let n = KdTreeNode{id: next_id, value: data[0].clone(), left: None, right: None, level: tree_depth+1, split_coordinate: 0 };
            return Some(Box::new(n));
        }

        sort_along_dimension(data, tree_depth % dimensionality);
        let median = data.len() / 2;
        let mut root = KdTreeNode{ id: next_id, value: data[median].clone(), left: None, right: None, level: tree_depth+1, split_coordinate: tree_depth % dimensionality };
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
/// use bioshell_clustering::kd_tree::{create_kd_tree, find_nearest};
/// use bioshell_clustering::{euclidean_distance_squared};
/// let mut points = [[0.1, 0.2], [0.2, 0.2], [1.1, 1.2], [2.2, 2.2]];
/// let root = create_kd_tree(&mut points, 2).unwrap();
/// let query = [0.3, 0.3];
/// let (d, e) = find_nearest(&root, &query, 2, euclidean_distance_squared);
/// assert!((d - 0.02).abs() < 0.000001);
/// ```
pub fn find_nearest<'a, T, F>(tree_root: &'a Box<KdTreeNode<T>>, query: &T, dimensionality: usize, distance:F) -> (f64, &'a T)
    where T: Index<usize, Output = f64>, T:Clone, F: Fn(&T, &T, usize) -> f64 {

    let mut min_dist: f64 = distance(tree_root.element(), query, dimensionality);
    let mut min_elem: &T = tree_root.element();
    let mut stack: Vec<&Box<KdTreeNode<T>>> = vec![];
    stack.push(tree_root);
    while !stack.is_empty() {
        let n = stack.pop().unwrap();
        let ne = n.element();
        let d = distance(ne, query, dimensionality);
        if d < min_dist { (min_dist, min_elem) = (d, ne); }

        let k = n.split_coordinate;
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
/// use bioshell_clustering::kd_tree::{create_kd_tree, find_within};
/// use bioshell_clustering::{euclidean_distance_squared};
/// let mut points = [[0.1, 0.2], [0.2, 0.2], [1.1, 1.2], [2.2, 2.2]];
/// let root = create_kd_tree(&mut points, 2).unwrap();
/// let query = [0.3, 0.3];
/// let neighbors = find_within(&root, &query, 2, 0.1, euclidean_distance_squared);
/// assert_eq!(neighbors.len(), 2);
/// ```
pub fn find_within<'a, T, F>(tree_root: &'a Box<KdTreeNode<T>>, query: &T, dimensionality: usize, radius: f64, distance:F) -> Vec<&'a T>
    where T: Index<usize, Output = f64>, T:Clone, F: Fn(&T, &T, usize) -> f64 {

    let mut ret: Vec<&T> = vec![];

    let mut stack: Vec<&Box<KdTreeNode<T>>> = vec![];
    stack.push(tree_root);
    while !stack.is_empty() {
        let n = stack.pop().unwrap();
        let ne = n.element();
        let d = distance(ne, query, dimensionality);
        if d < radius { ret.push(n.element()); }

        let k = n.split_coordinate;
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

/// Count points stored in a sub-tree rooted in a given node
///
/// # Examples
/// ```
/// # use rand::rngs::SmallRng;
/// # use rand::{SeedableRng, Rng};
/// # let mut rng = SmallRng::seed_from_u64(0);
/// use bioshell_clustering::kd_tree::{count, create_kd_tree};
/// let mut data = vec![vec![0.0]; 32];
/// for i in 0..data.len() { data[i][0] = rng.gen(); }
/// let root = create_kd_tree(&mut data, 1).unwrap();
/// assert_eq!(32, count(&root));
/// ```
pub fn count<T>(tree_node: &Box<KdTreeNode<T>>) -> usize where T: Index<usize, Output = f64> {

    fn add_count<T>(tree_node: &Box<KdTreeNode<T>>, cnt: usize) -> usize where T: Index<usize, Output = f64> {
        let mut sum = cnt;
        if let Some(left) = &tree_node.left { sum = add_count(&left, sum);}
        if let Some(right) = &tree_node.right { sum = add_count(&right, sum);}
        return sum + 1;
    }

    return add_count(tree_node,0);
}

fn sort_along_dimension<T>(data: &mut [T], idx:usize)  where T: Index<usize, Output = f64> {
    data.sort_by(|a,b|a[idx].partial_cmp(&b[idx]).unwrap());
}

#[cfg(test)]
mod test {
    use std::fmt::{Debug, Formatter};
    use rand::rngs::SmallRng;
    use rand::{Rng, SeedableRng};
    use crate::kd_tree::{create_kd_tree, KdTreeNode};

    #[test]
    fn test_kd_tree_construction() {
        const N: usize = 31;
        let mut rng = SmallRng::seed_from_u64(0);
        let mut data = vec![vec![0.0, 0.0]; N];
        for i in 0..data.len() {
            (data[i][0], data[i][1]) = (2.0 * rng.gen::<f64>() - 1.0, 2.0 * rng.gen::<f64>() - 1.0);
        }
        let mut root = create_kd_tree(&mut data.clone(), 2).unwrap();

        mark_partitions(&mut root,3,0);
        // println!("{:?}",&root);
        let mut min_v = vec![[2.0, 2.0]; 8];
        let mut max_v = vec![[-2.0, -2.0]; 8];
        for n in collect_nodes(&root) {
            min_v[n.level][0] = f64::min(min_v[n.level][0], n.value[0]);
            min_v[n.level][1] = f64::min(min_v[n.level][1], n.value[1]);
            max_v[n.level][0] = f64::max(max_v[n.level][0], n.value[0]);
            max_v[n.level][1] = f64::max(max_v[n.level][1], n.value[1]);
        }
        assert!(min_v[7][1] > max_v[6][1]);
        assert!(min_v[5][1] > max_v[4][1]);
        assert!(min_v[7][0] > max_v[5][0]);
        assert!(min_v[6][0] > max_v[4][0]);
    }

    fn collect_nodes<T>(root: &Box<KdTreeNode<T>>) -> Vec<&Box<KdTreeNode<T>>> {
        let mut ret: Vec<&Box<KdTreeNode<T>>> = vec![];
        ret.push(&root);

        fn collect_rec<'a, T>(root: &'a Box<KdTreeNode<T>>, store: &mut Vec<&'a Box<KdTreeNode<T>>>) {
            if let Some(left) = &root.left {
                store.push(left);
                collect_rec(left, store);
            }
            if let Some(right) = &root.right {
                store.push(right);
                collect_rec(right, store);
            }
        }
        collect_rec(root, &mut ret);
        return ret;
    }

    // KdTreeNode should be refactored into a tree data structure, then the following functions will
    // be much simplified.
    fn mark_partitions<T>(root: &mut Box<KdTreeNode<T>>, max_level: usize, last_mark: u32) {
        let mut mark = last_mark;
        if root.level <= max_level { mark = root.id; }

        if root.level > max_level { root.level = mark as usize; }
        if let Some(left) = &mut root.left {
            mark_partitions(left, max_level, mark);
        }
        if let Some(right) = &mut root.right {
            mark_partitions(right, max_level, mark);
        }
    }

    impl<T> Debug for KdTreeNode<T> where T: std::fmt::Debug{

        fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {

            fn fmt_rec<T>(node: &KdTreeNode<T>, out: & String) -> String where T: Debug {

                let mut o = format!("{}{}: {:?} {}\n", out, node.id, &node.value, node.level);
                if let Some(left) = &node.left {
                    o = fmt_rec(left, &o);
                }
                if let Some(right) = &node.right {
                    o = fmt_rec(right, &o);
                }
                return o;
            }

            let out = fmt_rec(&self,&String::default());
            write!(f, "{}", out)
        }
    }
}