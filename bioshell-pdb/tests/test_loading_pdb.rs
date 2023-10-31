#[cfg(test)]
mod tests {
    use std::io::BufReader;
    use std::string::String;
    use bioshell_pdb::{load_pdb_reader, PdbAtom};

    #[allow(non_upper_case_globals)]
    const pdb_2gb1:  &str = include_str!("./test_files/2gb1.pdb");

    #[test]
    fn test_2gb1_secondary_structure() {
        let strctr = load_pdb_reader(BufReader::new(pdb_2gb1.as_bytes())).unwrap();

        assert_eq!(strctr.count_atoms(), 855);
        assert_eq!(strctr.count_residues(), 56);
        assert_eq!(strctr.count_chains(), 1);


    }
}