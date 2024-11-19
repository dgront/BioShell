//! Generic binary tree struct and traversal operations
//!
//! BinaryTreeNode struct holds a value with the generic type T.
//! Option enum is used to represent its left and right child, which indeed are both optional.
//!
//! The following example builds a tree like this:
//! ```text
//!             1
//!          /    \
//!         2      3
//!        / \
//!       4   5
//! ```
//!
//! ```rust
//! use bioshell_datastructures::{BinaryTreeNode, depth_first_inorder};
//! // create a small tree
//! let root = BinaryTreeNode::new(1)
//!     .set_left(BinaryTreeNode::new(2)
//!             .set_left(BinaryTreeNode::new(4))
//!             .set_right(BinaryTreeNode::new(5))
//!     ).set_right(BinaryTreeNode::new(3));
//! assert!(root.left().is_some());
//! assert!(root.left().unwrap().value == 2);
//! let right = root.right().unwrap();
//! assert!(right.left().is_none());
//! assert!(right.right().is_none());
//!
//! // traverse the tree in-order to count its nodes
//! let mut node_cnt = 0;
//! depth_first_inorder(&root, &mut |_n| {node_cnt += 1});
//! assert_eq!(node_cnt, 5);
//! ```

/// Represents a node of a binary tree.
#[derive(Clone)]
pub struct BinaryTreeNode<T> {
    /// index of the node
    pub id: u32,
    /// the value stored on this node
    pub value: T,
    /// left sub-tree holds all points located to the left of the branching point
    pub left: Option<Box<BinaryTreeNode<T>>>,
    /// right sub-tree holds all points located to the right of the branching point
    pub right: Option<Box<BinaryTreeNode<T>>>,
}

impl<T> BinaryTreeNode<T> {

    /// Create a new node
    pub fn new(value: T) -> Self {
        BinaryTreeNode { id: 0, value, left: None, right: None }
    }

    /// Sets the left leaf of this node.
    ///
    /// Returns self for chained calls
    pub fn set_left<'a>(mut self, node: BinaryTreeNode<T>) -> Self {

        self.left = Some(Box::new(node));
        self
    }

    /// Sets the right leaf of this node.
    ///
    /// Returns self for chained calls
    pub fn set_right(mut self, node: BinaryTreeNode<T>) -> Self {
        self.right = Some(Box::new(node));
        self
    }

    /// Returns true if the node has a left child
    pub fn has_left( &self ) -> bool { self.left.is_some() }

    /// Returns true if the node has a right child
    pub fn has_right(&self) -> bool { self.right.is_some() }

    /// Returns true if the node has no children
    pub fn is_leaf(&self) -> bool { !self.has_left() && !self.has_right() }

    /// Borrows the left child of this node
    pub fn left(&self) -> Option<&BinaryTreeNode<T>> {
        self.left.as_ref().map(|n| n.as_ref())
    }

    /// Borrows the right child of this node
    pub fn right(&self) -> Option<&BinaryTreeNode<T>> {
        self.right.as_ref().map(|n| n.as_ref())
    }
}

/// Count all nodes of a given sub-tree, including the root node
///
/// # Examples
/// ```
/// use std::iter::zip;
/// use bioshell_datastructures::{BinaryTreeNode, collect_values, count_nodes};
/// let root = BinaryTreeNode::new(1)
///     .set_left(BinaryTreeNode::new(2)
///             .set_left(BinaryTreeNode::new(4))
///             .set_right(BinaryTreeNode::new(5))
///     ).set_right(BinaryTreeNode::new(3));
/// assert_eq!(count_nodes(&root), 5);
/// ```
pub fn count_nodes<T>(tree_node: &BinaryTreeNode<T>) -> usize {

    let mut ret: usize = 0;

    fn count_values_rec<'a, T>(tree_node: &'a BinaryTreeNode<T>, count: &mut usize) {
        if let Some(left) = &tree_node.left { count_values_rec(left, count);}
        *count += 1;
        if let Some(right) = &tree_node.right { count_values_rec(right, count);}
    }

    count_values_rec(tree_node, &mut ret);

    return ret;
}

/// Collect all data values stored in a sub-tree rooted in a given node
///
/// # Examples
/// ```
/// use std::iter::zip;
/// use bioshell_datastructures::{BinaryTreeNode, collect_values};
/// let root = BinaryTreeNode::new(1)
///     .set_left(BinaryTreeNode::new(2)
///             .set_left(BinaryTreeNode::new(4))
///             .set_right(BinaryTreeNode::new(5))
///     ).set_right(BinaryTreeNode::new(3));
/// let expected = vec![4, 2, 5, 1, 3];
/// for (a, b) in zip(&expected, collect_values(&root)) {
///     assert_eq!(a, b);
/// }
/// ```
pub fn collect_values<T>(tree_node: &BinaryTreeNode<T>) -> Vec<&T> {

    let mut ret: Vec<&T> = vec![];

    fn collect_values_rec<'a, T>(tree_node: &'a BinaryTreeNode<T>, leaf_elements: &mut Vec<&'a T>) {
        if let Some(left) = &tree_node.left { collect_values_rec(left, leaf_elements);}
        leaf_elements.push(&tree_node.value);
        if let Some(right) = &tree_node.right { collect_values_rec(right, leaf_elements);}
    }

    collect_values_rec(tree_node,&mut ret);

    return ret;
}

/// Collect data values stored in leaves of a sub-tree rooted in a given node
///
/// # Examples
/// ```
/// use std::iter::zip;
/// use bioshell_datastructures::{BinaryTreeNode, collect_leaf_values};
/// let root = BinaryTreeNode::new(1)
///     .set_left(BinaryTreeNode::new(2)
///             .set_left(BinaryTreeNode::new(4))
///             .set_right(BinaryTreeNode::new(5))
///     ).set_right(BinaryTreeNode::new(3));
/// let expected = vec![4, 5];
/// for (a, b) in zip(&expected, collect_leaf_values(&root)) {
///     assert_eq!(a, b);
/// }
/// ```
pub fn collect_leaf_values<T>(tree_node: &BinaryTreeNode<T>) -> Vec<&T> {

    let mut ret: Vec<&T> = vec![];

    fn collect_values_rec<'a, T>(tree_node: &'a BinaryTreeNode<T>, leaf_elements: &mut Vec<&'a T>) {
        if let Some(left) = &tree_node.left { collect_values_rec(left, leaf_elements);}
        if tree_node.is_leaf() { leaf_elements.push(&tree_node.value); }
        if let Some(right) = &tree_node.right { collect_values_rec(right, leaf_elements);}
    }

    collect_values_rec(tree_node,&mut ret);

    return ret;
}

/// Visits nodes of a binary tree in depth-first order, calls an action after seeing the left subtree.
///
/// # Examples
/// ```
/// // Count all nodes of a given tree
/// use bioshell_datastructures::{BinaryTreeNode, depth_first_inorder};
/// let root = BinaryTreeNode::new(1)
///     .set_left(BinaryTreeNode::new(2)
///             .set_left(BinaryTreeNode::new(4))
///             .set_right(BinaryTreeNode::new(5))
///     ).set_right(BinaryTreeNode::new(3));
/// let mut node_cnt = 0;
/// depth_first_inorder(&root, &mut |_n| {node_cnt += 1});
/// assert_eq!(node_cnt, 5);
/// ```
pub fn depth_first_inorder<T, F: FnMut(&BinaryTreeNode<T>)>(tree_node: &BinaryTreeNode<T>, action: &mut F) {

    fn inorder_rec<T, F: FnMut(&BinaryTreeNode<T>)>(tree_node: &BinaryTreeNode<T>, action: &mut F) {
        if let Some(left) = &tree_node.left { inorder_rec(left, action);}
        action(tree_node);
        if let Some(right) = &tree_node.right { inorder_rec(right, action);}
    }

    inorder_rec(tree_node, action);
}

pub fn depth_first_preorder<T, F: FnMut(&BinaryTreeNode<T>)>(tree_node: &BinaryTreeNode<T>, action: &mut F) {

    fn inorder_rec<T, F: FnMut(&BinaryTreeNode<T>)>(tree_node: &BinaryTreeNode<T>, action: &mut F) {
        action(tree_node);
        if let Some(left) = &tree_node.left { inorder_rec(left, action);}
        if let Some(right) = &tree_node.right { inorder_rec(right, action);}
    }

    inorder_rec(tree_node, action);
}
