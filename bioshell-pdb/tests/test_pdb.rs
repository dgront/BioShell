#[cfg(test)]
mod tests {
    use clap::Result;
    use bioshell_pdb::Pdb;

    #[test]
    fn test_pdb_from_file() -> Result<(), Box<dyn std::error::Error>> {
        let pdb_file_path = "tests/test_files/4hhb.pdb";
        let pdb = Pdb::from_file(pdb_file_path)?;

        assert_eq!(pdb.get_atoms_list().len(), 4384);

        let first_atom = &pdb.get_atoms_list()[0];
        assert_eq!(first_atom.get_atom_serial_no(), Some(0));
        assert_eq!(first_atom.get_atom_symbol(), "N");
        assert_eq!(first_atom.get_alt_loc_indicator(), " ");
        assert_eq!(first_atom.get_residue_name(), "VAL");
        assert_eq!(first_atom.get_chain_name(), "A");
        assert_eq!(first_atom.get_residue_no(), Some(1));
        assert_eq!(first_atom.get_insertion_code(), " ");
        assert_eq!(first_atom.get_coordinate().x, 11.54);
        assert_eq!(first_atom.get_coordinate().y, 11.88);
        assert_eq!(first_atom.get_coordinate().z, 7.95);
        assert_eq!(first_atom.get_occupancy(), Some(1.0));
        assert_eq!(first_atom.get_temperature_factor(), Some(0.0));
        assert_eq!(first_atom.get_atom_symbol(), "N");
        assert_eq!(first_atom.get_charge_of_the_atom(), " ");
        assert_eq!(first_atom.get_protein_name(), "4HHB");

        Ok(())
    }

    //#[test]
/*    fn test_pdb_write_csv() -> Result<(), Box<dyn std::error::Error>> {
        let pdb_file_path = "test_files/4HHB.pdb";
        let pdb = Pdb::from_file(pdb_file_path)?;
        let csv_file_path = "test_files/4HHB.csv";
        let expected_csv = include_str!("test_files/4HHB.csv");

        pdb.write_csv(csv_file_path)?;

        let csv_contents = std::fs::read_to_string(csv_file_path)?;
        assert_eq!(csv_contents, expected_csv);

        Ok(())
    }
    */
}