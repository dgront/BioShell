//! Generic binary tree struct and traversal operations
//!
//! BinaryTreeNode struct holds a value with the generic type T.
//! Option enum is used to represent its left and right child, which indeed are both optional.
//!
//! ```rust
//! use bioshell_datastructures::{BinaryTreeNode, depth_first_inorder};
//! // create a small tree
//! let root = BinaryTreeNode::new(1)
//!     .set_left(BinaryTreeNode::new(2)
//!             .set_left(BinaryTreeNode::new(4))
//!             .set_right(BinaryTreeNode::new(5))
//!     ).set_right(BinaryTreeNode::new(3));
//! // count its nodes
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
}

/// Count points stored in a sub-tree rooted in a given node
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
/// for (a, b) in zip(&expected, collect_values(&Box::new(root))) {
///     assert_eq!(a, b);
/// }
/// ```
pub fn collect_values<'a, T>(tree_node: &'a BinaryTreeNode<T>) -> Vec<&'a T> {

    let mut ret: Vec<&T> = vec![];

    fn collect_values_rec<'a, T>(tree_node: &'a BinaryTreeNode<T>, leaf_elements: &mut Vec<&'a T>) {
        if let Some(left) = &tree_node.left { collect_values_rec(left, leaf_elements);}
        leaf_elements.push(&tree_node.value);
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
pub fn depth_first_inorder<'a, T, F: FnMut(& BinaryTreeNode<T>)>(tree_node: &'a BinaryTreeNode<T>, action: &mut F) {

    fn inorder_rec<'a, T, F: FnMut(& BinaryTreeNode<T>)>(tree_node: &'a BinaryTreeNode<T>, action: &mut F) {
        if let Some(left) = &tree_node.left { inorder_rec(left, action);}
        action(tree_node);
        if let Some(right) = &tree_node.right { inorder_rec(right, action);}
    }

    inorder_rec(tree_node, action);
}
