#[cfg(test)]
mod tests {
    use bioshell_pdb::pdb_source::PdbSource;

    #[test]
    fn test_pdb_source() {
        let molecule_id = 1;
        let scientific_organism = "Homo sapiens";
        let common_organism = "human";
        let organism_tax_id = "9606";
        let strain = "ATCC";
        let pdb_source = PdbSource::new(molecule_id, scientific_organism, common_organism, organism_tax_id, strain);

        assert_eq!(pdb_source.get_molecule_id(), molecule_id);
        assert_eq!(pdb_source.get_scientific_organism(), scientific_organism);
        assert_eq!(pdb_source.get_common_organism(), common_organism);
        assert_eq!(pdb_source.get_organism_tax_id(), organism_tax_id);
        assert_eq!(pdb_source.get_strain(), strain);
    }
}