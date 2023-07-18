#[cfg(test)]
mod tests {
    use bioshell_pdb::pdb_helix_line_parser::PdbHelixParser;

    #[test]
    fn test_parse_line() {
        let line = "HELIX    1   A SER A   5  GLN A  11  1                                  7";
        let pdb_helix = PdbHelixParser::parse_line(line).unwrap();
        assert_eq!(pdb_helix.serial_number_1, 1);
        assert_eq!(pdb_helix.helix_id_2, "A");
        assert_eq!(pdb_helix.init_res_name_3, "SER");
        assert_eq!(pdb_helix.init_chain_id_4, 'A');
        assert_eq!(pdb_helix.init_seq_num_5, 5);
        assert_eq!(pdb_helix.init_insert_code_6, ' ');
        assert_eq!(pdb_helix.end_res_name_7, "GLN");
        assert_eq!(pdb_helix.end_chain_id_8, 'A');
        assert_eq!(pdb_helix.end_seq_num_9, 11);
        assert_eq!(pdb_helix.end_insert_code_10, ' ');
        assert_eq!(pdb_helix.helix_type_11, 1);
        assert_eq!(pdb_helix.helix_length_12, 7);
    }

    #[test]
    fn test_parse_line_invalid() {
        let line = "ATOM      1  N   GLY A   1       8.374  10.880  -6.778  1.00  0.00           N  ";
        let pdb_helix = PdbHelixParser::parse_line(line);
        assert!(pdb_helix.is_none());
    }
}