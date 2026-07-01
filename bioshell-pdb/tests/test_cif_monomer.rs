#[cfg(test)]
mod test_cif_monomer {

    #[test]
    fn cif_monomer_to_pdb() {
        let expected: [&str; 6] = [
            "ATOM      1  N   ALA A   1       2.281  26.213  12.804  1.00  0.00           N  ",
            "ATOM      2  CA  ALA A   1       1.169  26.942  13.411  1.00  0.00           C  ",
            "ATOM      3  C   ALA A   1       1.539  28.344  13.874  1.00  0.00           C  ",
            "ATOM      4  O   ALA A   1       2.709  28.647  14.114  1.00  0.00           O  ",
            "ATOM      5  CB  ALA A   1       0.601  26.143  14.574  1.00  0.00           C  ",
            "ATOM      6  OXT ALA A   1       0.523  29.194  13.997  1.00  0.00           O  ",
        ];

        let ala = include_str!("../tests/test_files/ala.cif");
        let reader = std::io::BufReader::new(ala.as_bytes());
        let strctr = bioshell_pdb::read_cif_monomers(reader, None).unwrap();
        let pdb_lines = strctr.atoms().iter().map(|a| a.to_string()).collect::<Vec<String>>();
        assert_eq!(pdb_lines.len(), 13);
        for i in 0..expected.len() {
            assert_eq!(pdb_lines[i], expected[i]);
        }
    }
}