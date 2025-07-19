#[cfg(test)]
mod tests {
    use bioshell_clustering::errors::ClusteringError;
    use bioshell_clustering::hierarchical::{balance_clustering_tree, retrieve_data, hierarchical_clustering, show_clustering_tree, DistanceMatrix, retrieve_clusters, medoid_by_min_max, retrieve_data_id};
    use bioshell_clustering::hierarchical::strategies::single_link;

    #[test]
    fn cluster_numbers() {

        let data: Vec<f32> = vec![6.0, 2.0, 2.5, 5.9];
        let distance_fn = |ia: usize, ib: usize| (data[ia] - data[ib]).abs();
        let mut clustering = hierarchical_clustering(data.len(), distance_fn, &single_link);
        assert_eq!(clustering.value.cluster_size, 4);
        let clustered_items = retrieve_data(&clustering, &data);
        assert_eq!(clustered_items, [6.0, 5.9, 2.0, 2.5]);

        balance_clustering_tree(&mut clustering, &distance_fn);
        let clustered_items = retrieve_data(&clustering, &data);
        assert_eq!(clustered_items, [6.0, 5.9, 2.5, 2.0]);

        println!("{}", show_clustering_tree(&clustering));
    }

    #[test]
    fn cluster_letters() {

        let data: Vec<char> = vec!['A', 'Z', 'D', 'Y', 'C', 'X', 'B', 'W', 'F', 'V'];
        let distance_fn = |ia: usize, ib: usize| (data[ia] as i16 - data[ib] as i16).abs() as f32;

        let mut clustering = hierarchical_clustering(data.len(), distance_fn, &single_link);
        assert_eq!(clustering.value.cluster_size, data.len());

        let clustered_items = retrieve_data(&clustering, &data);
        assert_eq!(clustered_items, ['A', 'B', 'D', 'C', 'F', 'Z', 'Y', 'X', 'W', 'V']);

        balance_clustering_tree(&mut clustering, &distance_fn);
        let clustered_items = retrieve_data(&clustering, &data);
        assert_eq!(clustered_items, ['A', 'B', 'C', 'D', 'F', 'V', 'W', 'X', 'Y', 'Z']);

        // println!("{}", show_clustering_tree(&clustering));
    }

    #[test]
    fn cluster_distance_matrix() -> Result<(), ClusteringError> {
        let dmatrix = DistanceMatrix::from_tsv("tests/test_files/d5.tsv")?;
        let distance_fn = |i: usize, j: usize| dmatrix.distance(i, j);
        let mut clustering = hierarchical_clustering(dmatrix.n_elements(), distance_fn, &single_link);
        let mut clusters = retrieve_clusters(&mut clustering, 0.11);
        clusters.sort_by(|a, b| a.value.cluster_size.cmp(&b.value.cluster_size));

        assert_eq!(clusters.len(), 4);
        let np = clusters.iter().map(|c|c.value.cluster_size).sum::<usize>();
        assert_eq!(np, 5);
        assert_eq!(clusters[3].value.cluster_size, 2);

        // --- the first 3 clusters have size 1 so the medoid must be equal to the sole cluster member
        for i in 0..3 {
            let medoid_idx = medoid_by_min_max(clusters[i], &distance_fn);
            let leaf_ids: Vec<usize> = retrieve_data_id(&clusters[i]);
            assert_eq!(medoid_idx, leaf_ids[0]);
        }

        Ok(())
    }
}