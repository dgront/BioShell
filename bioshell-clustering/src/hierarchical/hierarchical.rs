use std::collections::HashMap;
use std::fmt::{Display, Formatter};
use std::marker::PhantomData;
use std::ops::IndexMut;
use log::info;
use bioshell_datastructures::{BinaryTreeNode, collect_leaf_values, depth_first_inorder, depth_first_preorder};
use crate::hierarchical::{DistanceMatrix, single_link};

type  HierarchicalCluster = BinaryTreeNode<(usize, f32)>;

/// Defines a merging point in the hierarchical clustering process.
pub trait MergingPoint {
    /// size of a cluster, i.e. number of its nodes including itself
    fn size(&self) -> usize;

    /// distance at which this node has been created, i.e. the distance between the two leaves of this [`MergingPoint`]
    fn merging_distance(&self) -> f32;

    /// converts this [`MergingPoint`] to a string
    fn to_string(&self) -> String;
}

/// Extracts the [`MergingPoint`] data from the [`HierarchicalCluster`]
impl MergingPoint for HierarchicalCluster {

    fn size(&self) -> usize { self.value.0 }

    fn merging_distance(&self) -> f32 { self.value.1 }

    fn to_string(&self) -> String {
        let left_str = match self.left() {
            Some(l) => l.id.to_string(),
            None => "-".to_string()
        };
        let right_str = match self.right() {
            Some(l) => l.id.to_string(),
            None => "-".to_string()
        };
        format!("{} : {} {} @ {}", self.id, left_str, right_str, self.merging_distance())
    }
}

/// Provides hierarchical clustering calculations
pub struct HierarchicalClustering {
    n_points: usize,
    dmatrix: DistanceMatrix,
    clusters: HashMap<usize, HierarchicalCluster>,
    root: Option<HierarchicalCluster>,
}

impl HierarchicalClustering  {

    pub fn new<F: Fn(usize, usize) -> f32>(n_data: usize, distance: F) -> HierarchicalClustering {

        let dmatrix = DistanceMatrix::new(n_data, distance);
        let mut clusters = HashMap::new();
        for i in 0..n_data {
            let mut c = BinaryTreeNode::new((1, 0.0));
            c.id = i as u32;
            clusters.insert(i, c);
        }
        HierarchicalClustering{ n_points: n_data, dmatrix, clusters, root: None }
    }

    /// Calculates the hierarchical clustering
    pub fn cluster(&mut self) -> Option<&HierarchicalCluster> {

        let mut current_cluster_id = self.n_points;

        while self.clusters.len() > 1 {
            // --- find the two clusters with the smallest distance to be merged
            let (i, j) = self.dmatrix.closest_elements();   // it's guaranteed that i < j
            // --- find the two clusters themselves
            let ci = self.clusters.remove(&i).expect(format!("cluster {} not found", i).as_str());
            let cj = self.clusters.remove(&j).expect(format!("cluster {} not found", j).as_str());
            let id_i = ci.id;
            let id_j = cj.id;

            // --- merge the two clusters
            let merging_distance = self.dmatrix.matrix[i][j];
            let mut c = BinaryTreeNode::new((ci.size() + cj.size(), merging_distance));
            c = c.set_left(ci).set_right(cj);
            c.id = current_cluster_id as u32;
            info!("Merging clusters {} and {} into {} with distance {} at step {}",
                &id_i, &id_j, &c.id, merging_distance, current_cluster_id - self.n_points + 1);

            // --- the newly created cluster is inserted into the hashmap at "i"
            self.clusters.insert(i, c);
            // --- update the distance matrix for the newly created cluster
            self.dmatrix.update_distances(i, j, single_link, i);

            // --- the last cluster is moved to index "j" to keep the distance map consistent
            let last_cluster_idx = self.dmatrix.order - 1;
            if j < last_cluster_idx {
                let c = self.clusters.remove(&last_cluster_idx).expect(format!("cluster {} not found", last_cluster_idx).as_str());
                self.clusters.insert(j, c);
                // --- update the distance matrix for the moved cluster; this also lowers the order of the distance matrix by one
                self.dmatrix.replace_with_last(j);
            } else {
                self.dmatrix.order -= 1;
            }

            current_cluster_id += 1;
        }

        self.root = self.clusters.remove(&0);

        return self.root.as_ref();
    }

    /// borrows the root of the hierarchical clustering.
    ///
    /// The option is `None` if the clustering has not been calculated yet.
    pub fn root(&self) -> Option<&HierarchicalCluster> { self.root.as_ref() }

    pub fn clustered_data<T:Copy>(cluster: &HierarchicalCluster, all_data: &[T]) -> Vec<T> {

        let mut indexes: Vec<usize> = Vec::new();
        let mut collector = |n: &HierarchicalCluster| {
            if n.is_leaf() {
                indexes.push(n.id as usize);
            }
        };

        depth_first_preorder(cluster, &mut collector);
        let ret: Vec<T> = indexes.iter().map(|i| all_data[*i]).collect();

        return ret;
    }

    /// Balances the clustering tree to group most similar data items together.
    ///
    /// The method works by rotating recursively each subtree of the clustering tree starting from the bottom,
    /// trying to minimize the distance between the neighboring leaves.
    pub fn balance_clustering_tree<F: Fn(usize, usize) -> f32>(&mut self, distance: &F) {

        fn rotate_rec<F: Fn(usize, usize) -> f32>(tree_node: &mut HierarchicalCluster, distance: &F) {

            if let Some(left) = tree_node.left_mut() { rotate_rec(left, distance);}
            if let Some(right) = tree_node.right_mut() { rotate_rec(right, distance);}

            let rotate_what = HierarchicalClustering::if_rotate(tree_node, distance);
            if rotate_what.0 { tree_node.left_mut().unwrap().rotate(); }
            if rotate_what.1 { tree_node.right_mut().unwrap().rotate(); }
        }

        rotate_rec(self.root.as_mut().unwrap(), distance);
    }

    /// Finds the ID of the left-most leaf node.
    ///
    /// This method is used to balance the clustering tree
    fn get_leftmost_id(c: &HierarchicalCluster) -> usize {

        if c.has_left() {
            HierarchicalClustering::get_leftmost_id(c.left().unwrap())
        } else {
            c.id as usize
        }
    }

    /// Finds the ID of the right-most leaf node.
    ///
    /// This method is used to balance the clustering tree
    fn get_rightmost_id(c: &HierarchicalCluster) -> usize {

        if c.has_right() {
            HierarchicalClustering::get_rightmost_id(c.right().unwrap())
        } else {
            c.id as usize
        }
    }
    fn if_rotate<F: Fn(usize, usize) -> f32>(c: &mut HierarchicalCluster, distance: F)  -> (bool, bool) {

        if c.is_leaf() { return (false, false); }

        let left = c.left().unwrap();
        let right = c.right().unwrap();
        if right.is_leaf() && left.is_leaf() { return (false, false); }

        if right.is_leaf() {
            let left_right_idx = HierarchicalClustering::get_rightmost_id(left);
            let left_left_idx = HierarchicalClustering::get_leftmost_id(left);
            let right_idx = right.id as usize;
            if distance(right_idx, left_left_idx) < distance(right_idx, left_right_idx) {
                return (true, false);
            }
            return (false, false);
        }
        if left.is_leaf() {
            let right_right_idx = HierarchicalClustering::get_rightmost_id(right);
            let right_left_idx = HierarchicalClustering::get_leftmost_id(right);
            let left_idx = left.id as usize;
            if distance(left_idx, right_left_idx) > distance(left_idx, right_right_idx) {
                return (false, true);
            }
            return (false, false);
        }
        let right_right_idx = HierarchicalClustering::get_rightmost_id(right);
        let right_left_idx = HierarchicalClustering::get_leftmost_id(right);
        let left_right_idx = HierarchicalClustering::get_rightmost_id(left);
        let left_left_idx = HierarchicalClustering::get_leftmost_id(left);
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
}



impl Display for HierarchicalClustering {

    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {

        fn print_node_rec(n: &HierarchicalCluster, output_string: &mut String, indent: usize) {
            output_string.push_str(&"  ".repeat(indent));
            output_string.push_str(&n.to_string());
            output_string.push_str("\n");
            if n.has_left() {
                print_node_rec(n.left().unwrap(), output_string, indent + 1);
            }
            if n.has_right() {
                print_node_rec(n.right().unwrap(), output_string, indent + 1);
            }
        }
        let mut output_string: String = String::new();


        if let Some(root) = self.root() {
            print_node_rec(self.root().unwrap(), &mut output_string, 0);
        } else {
            output_string = "no clustering calculated yet".to_string();
        }

        write!(f, "{}", output_string)
    }
}
