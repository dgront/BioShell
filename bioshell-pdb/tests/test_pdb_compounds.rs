#[cfg(test)]
mod tests {
    use bioshell_pdb::pdb_compound::PdbCompound;

    #[test]
    fn test_pdb_compound() {
        let compound = PdbCompound::new(1, "my_molecule", "A");
        assert_eq!(compound.get_molecule_id(), 1);
        assert_eq!(compound.get_molecule_name(), "my_molecule");
        assert_eq!(compound.get_chain_id(), "A");
    }
}