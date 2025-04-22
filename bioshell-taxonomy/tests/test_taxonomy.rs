#[cfg(test)]
mod tests {
    use std::path::PathBuf;
    use bioshell_taxonomy::Taxonomy;

    const TEST_TAXDUMP_URL: &str = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz";
    const TEST_TAXDUMP_FILE: &str = "taxdump.tar.gz";

    #[test]
    fn test_file_loading() {
        // Download the file only if not already downloaded
        let path = PathBuf::from(TEST_TAXDUMP_FILE);
        // if !path.exists() {
        //     Taxonomy::download_taxdump_to_file(TEST_TAXDUMP_URL, &path)
        //         .expect("Failed to download taxdump");
        // }
        assert!(path.exists(), "Downloaded file does not exist");

        // Load taxonomy
        let taxonomy = Taxonomy::load_from_tar_gz(&path)
            .expect("Failed to load taxonomy");

        // Basic test: lookup Homo sapiens
        let human_taxid = taxonomy.taxid("Homo sapiens");
        assert!(human_taxid.is_some(), "Homo sapiens taxid not found");
        let human_taxid = human_taxid.unwrap();

        // Test species function
        let species_node = taxonomy.node(human_taxid);
        assert!(species_node.is_some(), "Species node for Homo sapiens not found");
        assert_eq!(species_node.unwrap().name, "Homo sapiens");

        // Test lineage function
        let lineage = taxonomy.lineage(human_taxid);
        assert!(!lineage.is_empty(), "Lineage for Homo sapiens should not be empty");
        assert_eq!(lineage.last().unwrap().name, "Homo sapiens");
    }


    #[test]
    fn test_basic_lookup() {
        let path = PathBuf::from(TEST_TAXDUMP_FILE);
        let taxonomy = Taxonomy::load_from_tar_gz(&path)
            .expect("Failed to load taxonomy");

        // Check if E. coli exists
        let ecoli_taxid = taxonomy.taxid("Escherichia coli");
        assert!(ecoli_taxid.is_some(), "E. coli taxid not found");

        let ecoli_node = taxonomy.node(ecoli_taxid.unwrap());
        assert!(ecoli_node.is_some(), "Species node for E. coli not found");
    }
}
