#[cfg(test)]
mod tests {
    use std::io::BufReader;
    use std::string::String;
    use bioshell_pdb::{Deposit};

    #[allow(non_upper_case_globals)]
    const pdb_2gb1:  &str = include_str!("./test_files/2gb1.pdb");

    #[test]
    fn load_2gb1_from_pdb() {
        let deposit = Deposit::from_pdb_reader(BufReader::new(pdb_2gb1.as_bytes())).unwrap();
        let strctr = deposit.structure().unwrap();

        assert_eq!(strctr.count_atoms(), 855);
        assert_eq!(strctr.count_residues(), 56);
        assert_eq!(strctr.count_chains(), 1);
        assert_eq!(strctr.count_models(), 1);
        let sec_str = strctr.secondary("A");
        assert_eq!(String::from("CEEEEEECCCCCCEEEEEECCHHHHHHHHHHHHHHHCCCCCEEEEECCCCEEEEEC"), sec_str.to_string());

        assert_eq!(deposit.title.unwrap().to_string(),
                   String::from("A NOVEL, HIGHLY STABLE FOLD OF THE IMMUNOGLOBULIN BINDING DOMAIN OF STREPTOCOCCAL PROTEIN G"))
    }
}