#[cfg(test)]
mod tests {
    use bioshell_clustering::hierarchical::{HierarchicalClustering, MergingPoint};
    use bioshell_datastructures::collect_leaf_values;

    #[test]
    fn cluster_numbers() {

        let data: Vec<f32> = vec![1.0, 2.0, 2.5, 4.9, 5.0, 8.3, 7.9, 8.0, 9.0, 9.4];
        let distance_fn = |ia: usize, ib: usize| (data[ia] - data[ib]).abs();

        let mut clustering = HierarchicalClustering::new(data.len(), distance_fn);
        assert!(clustering.root().is_none());

        let root = &clustering.cluster().unwrap();
        assert_eq!(root.size(), 10);

        let clustered_items = HierarchicalClustering::clustered_data(&root, &data);
        println!("clustered_items: {:?}", clustered_items);

        println!("{}", &clustering);
    }

    #[test]
    fn cluster_letters() {

        let data: Vec<u8> = vec![b'A', b'Z', b'B', b'Y', b'C', b'X', b'D', b'W', b'E', b'V', b'F'];
        let distance_fn = |ia: usize, ib: usize| (data[ia] as i16 - data[ib] as i16).abs() as f32;

        let mut clustering = HierarchicalClustering::new(data.len(), distance_fn);
        assert!(clustering.root().is_none());

        let root = &clustering.cluster().unwrap();
        assert_eq!(root.size(), 11);

        let clustered_items = HierarchicalClustering::clustered_data(&root, &data);
        println!("clustered_items: {:?}", clustered_items);

        println!("{}", &clustering);
    }
}