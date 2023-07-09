#[cfg(test)]
mod tests {
    use bioshell_pdb::pdb_sequence_of_residue::PdbSequenceOfResidue;

    #[test]
    fn test_pdb_sequence_of_residue() {
        let chain_id = "A";
        let sequence = vec!["ALA", "GLY", "LEU", "VAL", "PRO"];
        let pdb_sequence = PdbSequenceOfResidue::new(chain_id, sequence);

        assert_eq!(pdb_sequence.get_chain_id(), chain_id);
        assert_eq!(pdb_sequence.get_sequence(), &["ALA", "GLY", "LEU", "VAL", "PRO"]);
    }
}