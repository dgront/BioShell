#[cfg(test)]
mod tests {
    use bioshell_pdb::pdb_helix::PdbHelix;
    use bioshell_pdb::ResidueId;

    #[test]
    fn test_parse_line() {
        let line = "HELIX    1   A SER A    5  GLN A   11  1                                   7";
        let pdb_helix = PdbHelix::from_helix_line(line);
        assert_eq!(pdb_helix.ser_num, 1);
        assert_eq!(pdb_helix.helix_id, "A");
        assert_eq!(pdb_helix.init_res_name, "SER");
        assert_eq!(pdb_helix.init_chain_id, "A");
        assert_eq!(pdb_helix.init_seq_num, 5);
        assert_eq!(pdb_helix.init_i_code, ' ');
        assert_eq!(pdb_helix.end_res_name, "GLN");
        assert_eq!(pdb_helix.end_chain_id, "A");
        assert_eq!(pdb_helix.end_seq_num, 11);
        assert_eq!(pdb_helix.end_i_code, ' ');
        assert_eq!(pdb_helix.helix_class, 1);
        assert_eq!(pdb_helix.length, 7);

        assert_eq!(pdb_helix.init_res_id(), ResidueId::new("A", 5, ' '));
        assert_eq!(pdb_helix.end_res_id(), ResidueId::new("A", 11, ' '));
    }
}