use pyo3::prelude::*;
use pyo3::{Python, PyModule, PyClass, PyResult, pyclass, pymethods};

use bioshell_datastructures::KdTreeData<T>;

#[pyclass(name="BinaryTreeNode")]
pub struct PyKdTreeData<T> {
    data : KdTreeData<T>
}

#[pymethods]
impl<T> KdTreeData<T> {
    #[new]
    pub fn new(value: T) -> Self {
        KdTreeData{value: value, split_coordinate: 0}
    }
    pub fn level(&self) -> usize {
        self.level
    }
    pub fn split_coordinate(&self) -> usize {
        self.split_coordinate
    }
}


pub fn create_kd_tree<T>(data: &mut [T], dimensionality: usize) -> Option
 fn create_kd_tree_rec<T>(data: &mut [T], tree_depth: usize, dimensionality: usize, next_id: u32) -> Option
 fn sort_along_dimension<T>(data: &mut [T], dimension: usize)

pub fn find_nearest<T, F>(root: &BinaryTreeNode<KdTreeData<T>>, query: &[f64], dimensionality: usize, distance: F) -> (f64, &T)
    where T: Index<usize, Output = f64>, F: Fn(&[f64], &[f64]) -> f64 {
    let mut best_distance = std::f64::INFINITY;
    let mut best_point = &root.data.value;
    let mut current_node = root;
    let mut stack = Vec::new();
    loop {
        let current_distance = distance(&query, &current_node.data.value);
        if current_distance < best_distance {
            best_distance = current_distance;
            best_point = &current_node.data.value;
        }
        let split_coordinate = current_node.data.split_coordinate;
        let split_value = current_node.data.value[split_coordinate];
        let query_value = query[split_coordinate];
        let (near, far) = if query_value < split_value {
            (&current_node.left, &current_node.right)
        } else {
            (&current_node.right, &current_node.left)
        };
        if let Some(node) = near {
            stack.push(far);
            current_node = node;
        } else if let Some(node) = stack.pop() {
            current_node = node;
        } else {
            break;
        }
    }
    return (best_distance, best_point);
}

pub fn find_within <T, F>(root: &BinaryTreeNode<KdTreeData<T>>, query: &[f64], dimensionality: usize, distance: F, radius: f64) -> Vec<&T>
    where T: Index<usize, Output = f64>, F: Fn(&[f64], &[f64]) -> f64 {
    let mut result = Vec::new();
    let mut stack = Vec::new();
    stack.push(root);
    while let Some(current_node) = stack.pop() {
        let current_distance = distance(&query, &current_node.data.value);
        if current_distance < radius {
            result.push(&current_node.data.value);
        }
        let split_coordinate = current_node.data.split_coordinate;
        let split_value = current_node.data.value[split_coordinate];
        let query_value = query[split_coordinate];
        let (near, far) = if query_value < split_value {
            (&current_node.left, &current_node.right)
        } else {
            (&current_node.right, &current_node.left)
        };
        if let Some(node) = near {
            stack.push(node);
        }
        if let Some(node) = far {
            let split_distance = (query[split_coordinate] - current_node.data.value[split_coordinate]).abs();
            if split_distance < radius {
                stack.push(node);
            }
        }
    }
    return result;
}

pub fn sort_along_dimension<T>(data: &mut [T], dimension: usize) {
    data.sort_by(|a, b| a[dimension].partial_cmp(&b[dimension]).unwrap());
}


