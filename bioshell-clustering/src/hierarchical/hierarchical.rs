use std::collections::HashMap;
use log::info;
use bioshell_datastructures::{BinaryTreeNode, depth_first_preorder};
use crate::hierarchical::{DistanceMatrix};

/// Describes a single merging step in a hierarchical clustering process.
///
/// At each step of hierarchical clustering algorithm, two nearest clusters are found and merged into a new cluster.
/// The distance between the two merged clusters is termed the `merging_distance`, while `cluster_size`
/// is the number of data points in the new cluster.
///
pub struct HierarchicalCluster {
    /// The number of data points in this cluster
    pub cluster_size: usize,
    /// The distance at which this cluster was created
    pub merging_distance: f32,
}

type ClusteringTreeNode = BinaryTreeNode<HierarchicalCluster>;

pub fn hierarchical_clustering<F, M>(n_data: usize, distance_func: F,
            clustering_strategy: &M) -> ClusteringTreeNode
where
    F: Fn(usize, usize) -> f32,
    M: Fn(usize, usize, usize, f32, f32, f32) -> f32 {

    // --- create the distance matrix used by the clustering algorithm
    let mut dmatrix = DistanceMatrix::new(n_data, distance_func);

    // --- create the initial clusters
    let mut clusters = HashMap::new();
    for i in 0..n_data {
        let mut c = BinaryTreeNode::new(HierarchicalCluster { cluster_size: 1, merging_distance: 0.0 });
        c.id = i as u32;
        clusters.insert(i, c);
    }

    // --- run the clustering algorithm
    let mut current_cluster_id = n_data;
    while clusters.len() > 1 {
        // --- find the two clusters with the smallest distance to be merged
        let (i, j) = dmatrix.closest_elements();   // it's guaranteed that i < j
        // --- find the two clusters themselves
        let ci = clusters.remove(&i).expect(format!("cluster {} not found", i).as_str());
        let cj = clusters.remove(&j).expect(format!("cluster {} not found", j).as_str());
        let id_i = ci.id;
        let id_j = cj.id;

        // --- merge the two clusters
        let merging_distance = dmatrix.matrix[i][j];
        let new_merge = HierarchicalCluster{ cluster_size: ci.value.cluster_size + cj.value.cluster_size, merging_distance };
        let mut c = BinaryTreeNode::new(new_merge);
        c = c.set_left(ci).set_right(cj);
        c.id = current_cluster_id as u32;
        info!("Merging clusters {} and {} into {} with distance {} at step {}",
                &id_i, &id_j, &c.id, merging_distance, current_cluster_id - n_data + 1);

        // --- the newly created cluster is inserted into the hashmap at "i"
        clusters.insert(i, c);
        // --- update the distance matrix for the newly created cluster
        dmatrix.update_distances(i, j, clustering_strategy, i);

        // --- the last cluster is moved to index "j" to keep the distance map consistent
        let last_cluster_idx = dmatrix.order - 1;
        if j < last_cluster_idx {
            let c = clusters.remove(&last_cluster_idx).expect(format!("cluster {} not found", last_cluster_idx).as_str());
            clusters.insert(j, c);
            // --- update the distance matrix for the moved cluster; this also lowers the order of the distance matrix by one
            dmatrix.replace_with_last(j);
        } else {
            dmatrix.order -= 1;
        }

        current_cluster_id += 1;
    }

    return clusters.remove(&0).unwrap();
}

/// Balances the clustering tree to group most similar data items together.
///
/// The method works by rotating recursively each subtree of the clustering tree starting from the bottom,
/// trying to minimize the distance between the neighboring leaves.
pub fn balance_clustering_tree<F: Fn(usize, usize) -> f32>(root: &mut ClusteringTreeNode, distance: &F) {

    fn rotate_rec<F: Fn(usize, usize) -> f32>(tree_node: &mut ClusteringTreeNode, distance: &F) {

        if let Some(left) = tree_node.left_mut() { rotate_rec(left, distance);}
        if let Some(right) = tree_node.right_mut() { rotate_rec(right, distance);}

        let rotate_what = if_rotate(tree_node, distance);
        if rotate_what.0 { tree_node.left_mut().unwrap().rotate(); }
        if rotate_what.1 { tree_node.right_mut().unwrap().rotate(); }
    }

    rotate_rec(root, distance);
}

/// Copies the clustered data into a vector, according to the order defined by the clustering tree.
///
/// If you call [`balance_clustering_tree()`] before calling this method, the returned items will be grouped
/// together by their mutual distance.
pub fn clustered_data<T:Copy>(cluster: &ClusteringTreeNode, all_data: &[T]) -> Vec<T> {

    let mut indexes: Vec<usize> = Vec::new();
    let mut collector = |n: &ClusteringTreeNode| {
        if n.is_leaf() {
            indexes.push(n.id as usize);
        }
    };

    depth_first_preorder(cluster, &mut collector);
    let ret: Vec<T> = indexes.iter().map(|i| all_data[*i]).collect();

    return ret;
}

/// Finds the ID of the left-most leaf node.
///
/// This method is used to balance the clustering tree
fn get_leftmost_id(c: &ClusteringTreeNode) -> usize {

    if c.has_left() {
        get_leftmost_id(c.left().unwrap())
    } else {
        c.id as usize
    }
}

/// Finds the ID of the right-most leaf node.
///
/// This method is used to balance the clustering tree
fn get_rightmost_id(c: &ClusteringTreeNode) -> usize {

    if c.has_right() {
        get_rightmost_id(c.right().unwrap())
    } else {
        c.id as usize
    }
}
fn if_rotate<F: Fn(usize, usize) -> f32>(c: &mut ClusteringTreeNode, distance: F) -> (bool, bool) {

    if c.is_leaf() { return (false, false); }

    let left = c.left().unwrap();
    let right = c.right().unwrap();
    if right.is_leaf() && left.is_leaf() { return (false, false); }

    if right.is_leaf() {
        let left_right_idx = get_rightmost_id(left);
        let left_left_idx = get_leftmost_id(left);
        let right_idx = right.id as usize;
        if distance(right_idx, left_left_idx) < distance(right_idx, left_right_idx) {
            return (true, false);
        }
        return (false, false);
    }
    if left.is_leaf() {
        let right_right_idx = get_rightmost_id(right);
        let right_left_idx = get_leftmost_id(right);
        let left_idx = left.id as usize;
        if distance(left_idx, right_left_idx) > distance(left_idx, right_right_idx) {
            return (false, true);
        }
        return (false, false);
    }
    let right_right_idx = get_rightmost_id(right);
    let right_left_idx = get_leftmost_id(right);
    let left_right_idx = get_rightmost_id(left);
    let left_left_idx = get_leftmost_id(left);
    let distances = [
        distance(left_right_idx, right_left_idx),   // nothing to rotate
        distance(left_left_idx, right_left_idx),    // rotate left
        distance(left_right_idx, right_right_idx),  // rotate right
        distance(left_left_idx, right_right_idx),]; // rotate both
    let index_of_min: usize = distances
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.total_cmp(b))
        .map(|(index, _)| index).unwrap();
    match index_of_min {
        1 => return (true, false),
        2 => return (false, true),
        3 => return (true, true),
        _ => return (false, false)
    }
}

pub fn show_clustering_tree(root: &ClusteringTreeNode) -> String {
    fn print_node_rec(n: &ClusteringTreeNode, output_string: &mut String, indent: usize) {
        output_string.push_str(&"  ".repeat(indent));
        output_string.push_str(&node_to_string(n));
        output_string.push_str("\n");
        if n.has_left() {
            print_node_rec(n.left().unwrap(), output_string, indent + 1);
        }
        if n.has_right() {
            print_node_rec(n.right().unwrap(), output_string, indent + 1);
        }
    }
    let mut output_string: String = String::new();

    print_node_rec(root, &mut output_string, 0);

    return output_string;
}

fn node_to_string(n: &ClusteringTreeNode) -> String {
    let left_str = match n.left() {
        Some(l) => l.id.to_string(),
        None => "-".to_string()
    };
    let right_str = match n.right() {
        Some(l) => l.id.to_string(),
        None => "-".to_string()
    };
    format!("{} : {} {} @ {}", n.id, left_str, right_str, n.value.merging_distance)
}
