#[cfg(test)]
mod tests {
    use bioshell_clustering::hierarchical::{balance_clustering_tree, clustered_data, hierarchical_clustering, show_clustering_tree};
    use bioshell_clustering::hierarchical::strategies::single_link;

    #[test]
    fn cluster_numbers() {

        let data: Vec<f32> = vec![6.0, 2.0, 2.5, 5.9];
        let distance_fn = |ia: usize, ib: usize| (data[ia] - data[ib]).abs();
        let mut clustering = hierarchical_clustering(data.len(), distance_fn, single_link);
        assert_eq!(clustering.value.cluster_size, 4);
        let clustered_items = clustered_data(&clustering, &data);
        assert_eq!(clustered_items, [6.0, 5.9, 2.0, 2.5]);

        balance_clustering_tree(&mut clustering, &distance_fn);
        let clustered_items = clustered_data(&clustering, &data);
        assert_eq!(clustered_items, [6.0, 5.9, 2.5, 2.0]);

        println!("{}", show_clustering_tree(&clustering));
    }

    #[test]
    fn cluster_letters() {

        let data: Vec<char> = vec!['A', 'Z', 'D', 'Y', 'C', 'X', 'B', 'W', 'F', 'V'];
        let distance_fn = |ia: usize, ib: usize| (data[ia] as i16 - data[ib] as i16).abs() as f32;

        let mut clustering = hierarchical_clustering(data.len(), distance_fn, single_link);
        assert_eq!(clustering.value.cluster_size, data.len());

        let clustered_items = clustered_data(&clustering, &data);
        assert_eq!(clustered_items, ['A', 'B', 'D', 'C', 'F', 'Z', 'Y', 'X', 'W', 'V']);

        balance_clustering_tree(&mut clustering, &distance_fn);
        let clustered_items = clustered_data(&clustering, &data);
        assert_eq!(clustered_items, ['A', 'B', 'C', 'D', 'F', 'V', 'W', 'X', 'Y', 'Z']);

        // println!("{}", show_clustering_tree(&clustering));
    }
}