use std::iter::zip;
use bioshell_datastructures::{BinaryTreeNode, collect_values, depth_first_inorder};

#[cfg(feature = "vec3")]
use bioshell_numerical::Vec3;

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


#[cfg(test)]
mod test_kdtree {
    use bioshell_datastructures::{depth_first_inorder, BinaryTreeNode, collect_values};
    use bioshell_datastructures::kd_tree::{euclidean_distance_squared, create_kd_tree, find_nearest, find_within, KdTreeData};

    use rand::rngs::SmallRng;
    use rand::SeedableRng;
    use rand::Rng;

    #[test]
    fn test_kd_tree() {
        let mut data = vec![vec![0.0, 0.5], vec![0.7, 0.5], vec![0.8, 1.1]];
        create_kd_tree(&mut data, 2);
    }


    #[test]
    fn test_k1_tree() {
        const N: usize = 128;
        let mut rng = SmallRng::seed_from_u64(0);
        let mut data = vec![vec![0.0]; N];
        for i in 0..data.len() { data[i][0] = rng.gen(); }

        let query = vec![0.5];
        // --- find nearest with brute force
        let (mut min_d, mut min_e) = (euclidean_distance_squared(&query, &data[0], 1), &data[0]);
        for e in data.iter() {
            let d = euclidean_distance_squared(&query, e, 1);
            if d < min_d { (min_d, min_e) = (d, e); }
        }

        let root = create_kd_tree(&mut data.clone(), 1).unwrap();
        let mut node_cnt = 0;
        depth_first_inorder(&root, &mut |_n| {node_cnt += 1});
        assert_eq!(N, node_cnt);

        // --- find nearest with KdTree
        let (d, e) = find_nearest(&root, &query, 1, euclidean_distance_squared);
        // --- check if we have the same element
        assert!((d - min_d).abs() < 0.000001);
        assert!((e[0] - min_e[0]).abs() < 0.000001);

        let nbors = find_within(&root, &query, 1, 0.1 * 0.1, euclidean_distance_squared);
        for e in nbors { assert!(euclidean_distance_squared(e, &query, 1) <= 0.1 * 0.1) }
    }

    /// Use ``cargo test --features "vec3"`` command to run these two tests:
    #[test]
    #[allow(non_snake_case)]
    #[cfg(feature = "vec3")]
    fn test_kd_tree_Vec3_2D() {
        const N: usize = 7;
        let mut data = vec![Vec3::new(0.0, 0.0, 0.0); N * N];
        for i in 0..N * N {
            data[i].x = (i % N) as f64 * 0.1 + 0.1;
            data[i].y = (i / N) as f64 * 0.1 + 0.1;
        }
        let root = create_kd_tree(&mut data.clone(), 2).unwrap();
        let neighbors = find_within(&root, &Vec3::new(0.3, 0.3, 0.0),
                                    2, 0.021, euclidean_distance_squared);
        assert_eq!(neighbors.len(), 9);
    }

    #[test]
    #[allow(non_snake_case)]
    #[cfg(feature = "vec3")]
    fn test_kd_tree_Vec3_3D() {
        const N: usize = 10;
        let mut data = vec![Vec3::new(0.0, 0.0, 0.0); N * N * N];
        for k in 0..N {
            for j in 0..N {
                for i in 0..N {
                    let idx = N * N * k + N * j + i;
                    data[idx].x = i as f64 * 0.1 + 0.1;
                    data[idx].y = j as f64 * 0.1 + 0.1;
                    data[idx].z = k as f64 * 0.1 + 0.1;
                }
            }
        }
        let root = create_kd_tree(&mut data.clone(), 3).unwrap();
        let neighbors = find_within(&root, &Vec3::new(0.3, 0.3, 0.3),
                                    3, 0.031, euclidean_distance_squared);

        assert_eq!(neighbors.len(), 27);
    }

    fn mark_partitions<T>(root: &mut Box<BinaryTreeNode<KdTreeData<T>>>, max_level: usize, last_mark: u32) {
        let mut mark = last_mark;
        if root.value.level <= max_level { mark = root.id; }

        if root.value.level > max_level { root.value.level = mark as usize; }
        if let Some(left) = &mut root.left {
            mark_partitions(left, max_level, mark);
        }
        if let Some(right) = &mut root.right {
            mark_partitions(right, max_level, mark);
        }
    }

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
        for n in collect_values(&root) {
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
}
