use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

#[pyclass]
#[derive(Clone)]
pub struct BinaryTreeNode<T> {
    #[pyo3(get, set)]
    pub id: u32,
    #[pyo3(get, set)]
    pub value: T,
    pub left: Option<Box<BinaryTreeNode<T>>>,
    pub right: Option<Box<BinaryTreeNode<T>>>,
}

#[pymethods]
impl<T> BinaryTreeNode<T> {
    #[new]
    pub fn new(value: T) -> Self {
        BinaryTreeNode {
            id: 0,
            value,
            left: None,
            right: None,
        }
    }

    pub fn set_left(&mut self, node: BinaryTreeNode<T>) -> PyResult<()> {
        self.left = Some(Box::new(node));
        Ok(())
    }

    pub fn set_right(&mut self, node: BinaryTreeNode<T>) -> PyResult<()> {
        self.right = Some(Box::new(node));
        Ok(())
    }
}

#[pyfunction]
pub fn collect_values<'a, T>(tree_node: &'a BinaryTreeNode<T>) -> Vec<&'a T> {
    let mut ret: Vec<&T> = vec![];

    fn collect_values_rec<'a, T>(tree_node: &'a BinaryTreeNode<T>, leaf_elements: &mut Vec<&'a T>) {
        if let Some(left) = &tree_node.left { collect_values_rec(left, leaf_elements);}
        leaf_elements.push(&tree_node.value);
        if let Some(right) = &tree_node.right { collect_values_rec(right, leaf_elements);}
    }

    collect_values_rec(tree_node, &mut ret);

    ret
}

#[pyfunction]
pub fn depth_first_inorder<'a, T, F>(tree_node: &'a BinaryTreeNode<T>, mut action: F) 
where 
    F: FnMut(&BinaryTreeNode<T>) 
{
    fn inorder_rec<'a, T, F>(tree_node: &'a BinaryTreeNode<T>, action: &mut F)
    where
        F: FnMut(&BinaryTreeNode<T>)
    {
        if let Some(left) = &tree_node.left { inorder_rec(left, action);}
        action(tree_node);
        if let Some(right) = &tree_node.right { inorder_rec(right, action);}
    }

    inorder_rec(tree_node, &mut action);
}

#[pymodule]
fn rust_binary_tree(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<BinaryTreeNode<i32>>()?;
    m.add_function(wrap_pyfunction!(collect_values, m)?)?;
    m.add_function(wrap_pyfunction!(depth_first_inorder, m)?)?;
    Ok(())
}