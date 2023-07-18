#[cfg(test)]
mod tests {
    use bioshell_pdb::pdb_sheet_line_parser::PdbSheetParser;

    #[test]
    fn test_parse_line() {
        let line = "SHEET    1   A 2 LEU A   3  ARG A   6 -1  N  LEU A   3   O  ARG A   6           ";
        let pdb_sheet = PdbSheetParser::parse_line(line).unwrap();
        assert_eq!(pdb_sheet.strand_1, 1);
        assert_eq!(pdb_sheet.sheet_id_2, "A");
        assert_eq!(pdb_sheet.num_strands_3, 2);
        assert_eq!(pdb_sheet.init_res_name_4, "LEU");
        assert_eq!(pdb_sheet.init_chain_id_5, 'A');
        assert_eq!(pdb_sheet.init_seq_num_6, 3);
        assert_eq!(pdb_sheet.init_insert_code_7, ' ');
        assert_eq!(pdb_sheet.end_res_name_8, "ARG");
        assert_eq!(pdb_sheet.end_chain_id_9, 'A');
        assert_eq!(pdb_sheet.end_seq_num_10, 6);
        assert_eq!(pdb_sheet.end_insert_code_11, ' ');
        assert_eq!(pdb_sheet.sense_12, -1);
    }

    #[test]
    fn test_parse_line_invalid() {
        let line = "ATOM      1  N   GLY A   1       8.374  10.880  -6.778  1.00  0.00           N  ";
        let pdb_sheet = PdbSheetParser::parse_line(line);
        assert!(pdb_sheet.is_none());
    }
}