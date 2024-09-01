use pyo3::prelude::*;
use pyo3::{Python, PyModule, PyClass, PyResult, pyclass, pymethods};

use bioshell_datastructures::BinaryTreeNode;

#[pyclass(name="BinaryTreeNode")]
pub struct PyBinaryTreeNode<T> {
    data : BinaryTreeNode<T>
}

#[pymethods]
impl<T> PyBinaryTreeNode<T> {
    #[new]
    pub fn new(value: T) -> Self {
        PyBinaryTreeNode{data: BinaryTreeNode::new(value)}
    }

    pub fn set_left(&mut self, node: PyBinaryTreeNode<T>){
        self.data.set_left(node.data);
    }

    pub fn set_right(&mut self, node: PyBinaryTreeNode<T>){
        self.data.set_right(node.data);
    }

    pub fn collect_values(&self) -> Vec<&T>{
        self.data.collect_values()
    }

    pub fn depth_first_inorder(&self) -> Vec<&T>{
        self.data.depth_first_inorder()
    }
}