#[cfg(test)]
mod tests {
    use std::error::Error;
    use std::path::PathBuf;
    use bioshell_taxonomy::{Rank, Taxonomy, TaxonomyMatcher};

    const TEST_TAXDUMP_FILE: &str = "./tests/test_files/test_taxdump.tar.gz";

    #[test]
    fn test_file_loading() {
        let path = PathBuf::from(TEST_TAXDUMP_FILE);
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
        assert_eq!(lineage.len(), 32, "Lineage for Homo sapiens should have 32 nodes");
        assert_eq!(lineage.last().unwrap().name, "Homo sapiens");
        lineage.iter().any(|node| node.name == "Vertebrate");
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

    #[test]
    fn test_taxonomy_matcher() -> Result<(), Box<dyn Error>> {
        let path = PathBuf::from(TEST_TAXDUMP_FILE);
        let taxonomy = Taxonomy::load_from_tar_gz(&path)
            .expect("Failed to load taxonomy");
        let matcher = TaxonomyMatcher::from_taxonomy(taxonomy)?;

        for name in ["Mus musculus", "mouse", "house mouse"] {
            let taxid = matcher.find(name);
            assert!(taxid.is_some(), "Taxid for {} not found", name);
        }
        Ok(())
    }

    #[test]
    fn list_kingdoms() -> Result<(), Box<dyn Error>> {
        // let path = PathBuf::from("taxdump.tar.gz");
        let path = PathBuf::from(TEST_TAXDUMP_FILE);
        let taxonomy = Taxonomy::load_from_tar_gz(&path)
            .expect("Failed to load taxonomy from taxdump.tar.gz");

        let mut kingdoms = vec![];
        for n in taxonomy.nodes().filter(|ni|ni.rank==Rank::Kingdom) {
            kingdoms.push(n.tax_id);
            // let name = &n.name;
            // print!("{} {} : ",n.tax_id, &name);
            // for synonym in taxonomy.names(n.tax_id) {
            //     if synonym != name { print!("{}; ", synonym); }
            // }
            // println!()
        }
        assert_eq!(kingdoms.len(), 2);
        Ok(())
    }

    #[test]
    fn test_serialize() -> Result<(), Box<dyn Error>> {
        let path = PathBuf::from(TEST_TAXDUMP_FILE);
        let taxonomy = Taxonomy::load_from_tar_gz(&path)
            .expect("Failed to load taxonomy from taxdump.tar.gz");

        let ecoli_taxid = taxonomy.taxid("Escherichia coli").expect("E. coli taxid not found");
        let node = taxonomy.node(ecoli_taxid);
        let json = serde_json::to_string_pretty(&node).expect("Can't serialize a Node struct!");

        let expected = r#"{
  "tax_id": 562,
  "parent_tax_id": 561,
  "rank": "Species",
  "name": "Escherichia coli"
}"#;
        assert_eq!(&json, expected);
        Ok(())
    }
}
