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

        assert_eq!(pdb.get_atoms_list().len(), 3801);
        //"ATOM      1  N   VAL A   1      19.323  29.727  42.781  1.00 49.05           N  " //4hhb
        //ATOM      1  N   GLU A   5     -15.953  21.156  16.122  1.00 30.79           N  "//16pk
        let first_atom = &pdb.get_atoms_list()[0];
        assert_eq!(first_atom.atom_serial_no_1, Some(1));
        assert_eq!(first_atom.atom_symbol_2, "N");
        assert_eq!(first_atom.atom_position_3, "");
        assert_eq!(first_atom.atom_no_in_the_branch_4, Some(0));
        assert_eq!(first_atom.connected_to_atom_no_5, Some(0));
        assert_eq!(first_atom.alt_loc_indicator_6, "");
        assert_eq!(first_atom.residue_name_7, "GLU");
        assert_eq!(first_atom.chain_name_8, "A");
        assert_eq!(first_atom.residue_no_9, Some(5));
        assert_eq!(first_atom.insertion_code_10, "");
        assert_eq!(first_atom.coordinate_11.x, -15.953);
        assert_eq!(first_atom.coordinate_11.y, 21.156);
        assert_eq!(first_atom.coordinate_11.z, 16.122);
        assert_eq!(first_atom.occupancy_12, Some(1.0));
        assert_eq!(first_atom.temperature_factor_13, Some(30.79));
        assert_eq!(first_atom.segment_identifier_14, "");
        assert_eq!(first_atom.segment_identifier_symbol_15, "N");
        assert_eq!(first_atom.charge_of_the_atom_16, "");
        assert_eq!(first_atom.protein_name, "16PK");

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