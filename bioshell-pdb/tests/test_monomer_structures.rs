#[cfg(test)]
mod tests {
    use bioshell_pdb::monomers::MonomerManager;

    #[test]
    fn test_monomer_manager() {

        let manager = MonomerManager::get();
        assert!(manager.by_code3("ALA").is_some());
        assert_eq!(manager.count(), 29);
    }

    #[test]
    fn test_ala_structure() {

        let manager = MonomerManager::get();
        assert!(manager.by_code3("ALA").is_some());
        let ala = manager.by_code3("ALA").unwrap();
        assert_eq!(ala.count_all_atoms(), 13);
        assert_eq!(ala.count_residue_atoms(), 10);
        assert_eq!(ala.count_residue_heavy(), 5);
    }
}