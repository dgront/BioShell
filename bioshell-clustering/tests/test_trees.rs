use std::iter::zip;
use bioshell_datastructures::{BinaryTreeNode, collect_values, depth_first_inorder};

#[test]
fn create_tree() {

    let mut tree: BinaryTreeNode<usize> = BinaryTreeNode::new(1);
    assert_eq!(tree.value, 1);
    tree = tree.set_left(BinaryTreeNode::new(2));
    if let Some(left) = &tree.left { assert_eq!(left.value, 2); }

    tree = tree.set_right(BinaryTreeNode::new(3));
    if let Some(right) = &tree.right { assert_eq!(right.value, 3); }

    let mut node_cnt = 0;
    depth_first_inorder(&tree, &mut |_n| { node_cnt += 1 });
    assert_eq!(node_cnt, 3);
}

#[test]
fn check_collected_values() {
    let root = BinaryTreeNode::new(1)
        .set_left(BinaryTreeNode::new(2)
            .set_left(BinaryTreeNode::new(4))
            .set_right(BinaryTreeNode::new(5))
        ).set_right(BinaryTreeNode::new(3));
    let expected = vec![4, 2, 5, 1, 3];
    for (a, b) in zip(&expected, collect_values(&Box::new(root))) {
        assert_eq!(a, b);
    }
}