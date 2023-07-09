#[cfg(test)]
mod tests {
    use clap::Result;
    use bioshell_pdb::Pdb;

    #[test]
    fn test_pdb_from_file() -> Result<(), Box<dyn std::error::Error>> {
        let pdb_file_path = "test_files/4HHB.pdb";
        let pdb = Pdb::from_file(pdb_file_path)?;

        assert_eq!(pdb.atoms_list.len(), 4770);

        let first_atom = &pdb.atoms_list[0];
        assert_eq!(first_atom.atom_serial_number, 1);
        assert_eq!(first_atom.atom_name, "N");
        assert_eq!(first_atom.alternate_location_indicator, ' ');
        assert_eq!(first_atom.residue_name, "VAL");
        assert_eq!(first_atom.chain_identifier, 'A');
        assert_eq!(first_atom.residue_sequence_number, 1);
        assert_eq!(first_atom.insertion_code, ' ');
        assert_eq!(first_atom.coordinate_x, 11.54);
        assert_eq!(first_atom.coordinate_y, 11.88);
        assert_eq!(first_atom.coordinate_z, 7.95);
        assert_eq!(first_atom.occupancy, 1.0);
        assert_eq!(first_atom.temperature_factor, 0.0);
        assert_eq!(first_atom.element_symbol, "N");
        assert_eq!(first_atom.charge, " ");
        assert_eq!(first_atom.protein_name, "4HHB");

        Ok(())
    }

    #[test]
    fn test_pdb_write_csv() -> Result<(), Box<dyn std::error::Error>> {
        let pdb_file_path = "test_files/4HHB.pdb";
        let pdb = Pdb::from_file(pdb_file_path)?;
        let csv_file_path = "test_files/4HHB.csv";
        let expected_csv = include_str!("../test_files/4HHB.csv");

        pdb.write_csv(csv_file_path)?;

        let csv_contents = std::fs::read_to_string(csv_file_path)?;
        assert_eq!(csv_contents, expected_csv);

        Ok(())
    }
}