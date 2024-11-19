#[cfg(test)]
mod tests {
    use bioshell_clustering::hierarchical::{HierarchicalClustering, MergingPoint};

    #[test]
    fn cluster_numbers() {

        let data: Vec<f32> = vec![6.0, 2.0, 2.5, 5.9];
        let distance_fn = |ia: usize, ib: usize| (data[ia] - data[ib]).abs();

        let mut clustering = HierarchicalClustering::new(data.len(), distance_fn);
        assert!(clustering.root().is_none());

        clustering.cluster();
        assert_eq!(clustering.root().unwrap().size(), 4);
        let clustered_items = HierarchicalClustering::clustered_data(clustering.root().unwrap(), &data);
        assert_eq!(clustered_items, [6.0, 5.9, 2.0, 2.5]);

        clustering.balance_clustering_tree(&distance_fn);
        let clustered_items = HierarchicalClustering::clustered_data(clustering.root().unwrap(), &data);
        assert_eq!(clustered_items, [6.0, 5.9, 2.5, 2.0]);
    }

    #[test]
    fn cluster_letters() {

        let data: Vec<char> = vec!['A', 'Z', 'D', 'Y', 'C', 'X', 'B', 'W', 'F', 'V'];
        let distance_fn = |ia: usize, ib: usize| (data[ia] as i16 - data[ib] as i16).abs() as f32;

        let mut clustering = HierarchicalClustering::new(data.len(), distance_fn);
        assert!(clustering.root().is_none());

        let root = clustering.cluster().unwrap();
        assert_eq!(root.size(), data.len());

        let clustered_items = HierarchicalClustering::clustered_data(clustering.root().unwrap(), &data);
        assert_eq!(clustered_items, ['A', 'B', 'D', 'C', 'F', 'Z', 'Y', 'X', 'W', 'V']);

        clustering.balance_clustering_tree(&distance_fn);
        let clustered_items = HierarchicalClustering::clustered_data(clustering.root().unwrap(), &data);
        assert_eq!(clustered_items, ['A', 'B', 'C', 'D', 'F', 'V', 'W', 'X', 'Y', 'Z']);

        // println!("{}", &clustering);
    }
}