#[cfg(test)]
mod tests {
    use std::env;
    use clap::Result;
    use bioshell_pdb::Pdb;

    #[test]
    fn test_pdb_from_file() -> Result<(), Box<dyn std::error::Error>> {
        if let Ok(current_dir) = env::current_dir() {
            println!("Current working directory: {}", current_dir.display());
        }
        let pdb_file_path = "./tests/test_files/16pk.pdb";
        let pdb = Pdb::from_file(pdb_file_path)?;

        assert_eq!(pdb.get_atoms_list().len(), 3781);
        //"ATOM      1  N   VAL A   1      19.323  29.727  42.781  1.00 49.05           N  "
        let first_atom = &pdb.get_atoms_list()[0];
        assert_eq!(first_atom.atom_serial_no_1, Some(0));
        assert_eq!(first_atom.atom_symbol_2, "N");
        assert_eq!(first_atom.alt_loc_indicator_6, " ");
        assert_eq!(first_atom.residue_name_7, "VAL");
        assert_eq!(first_atom.chain_name_8, "A");
        assert_eq!(first_atom.residue_no_9, Some(1));
        assert_eq!(first_atom.insertion_code_10, " ");
        assert_eq!(first_atom.coordinate_11.x, 11.54);
        assert_eq!(first_atom.coordinate_11.y, 11.88);
        assert_eq!(first_atom.coordinate_11.z, 7.95);
        assert_eq!(first_atom.occupancy_12, Some(1.0));
        assert_eq!(first_atom.temperature_factor_13, Some(0.0));
        assert_eq!(first_atom.atom_symbol_2, "N");
        assert_eq!(first_atom.charge_of_the_atom_16, " ");
        assert_eq!(first_atom.protein_name, "4HHB");

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