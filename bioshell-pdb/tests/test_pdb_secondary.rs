#[cfg(test)]
mod tests {
    use std::io::BufReader;
    use bioshell_pdb::{Deposit, SecondaryStructure};

    #[allow(non_upper_case_globals)]
    const cif_2gb1:  &str = include_str!("./test_files/2gb1.cif");
    #[allow(non_upper_case_globals)]
    const cif_2fdo: &str = include_str!("./test_files/2fdo.cif");
    #[allow(non_upper_case_globals)]
    const cif_4esa: &str = include_str!("./test_files/4esa.cif");

    #[test]
    fn secondary_from_cif() {
        let test_cases = vec![(cif_2gb1, "CEEEEEECCCCCCEEEEEECCHHHHHHHHHHHHHHHCCCCCEEEEECCCCEEEEEC"),
            (cif_2fdo, "CCCCCCEEEHHHHHHHHHHHHHCCCEEEEEECCEEEEEEEEECCCCCCEEEEEEEEEEECCCCCCHHHHHHHHCEEEEEEEECHHHHHHHHHH"),
            (cif_4esa, "CCCCHHHHHHHHHHHHHHHHCHHHHHHHHHHHHHHHHHHHHHHHHHCCCCCCCCHHHHHHHHHHHHHHHHHHHHHCCHHHHHHHHHHHHHHHCCCCHHHHHHHHHHHHHHHHHHHHCCCCHHHHHHHHHHHHHHHHHHHHCCC")];

        for (input, expctd) in test_cases {
            let reader = BufReader::new(input.as_bytes());
            let deposit = Deposit::from_cif_reader(reader);
            assert!(deposit.is_ok());
            let strctr = deposit.unwrap().structure().unwrap();

            assert_eq!(strctr.secondary("A").to_string(), expctd);
        }
    }

    #[test]
    fn create_secondary_structure() {
        let sec_str = SecondaryStructure::new(10);
        assert_eq!(sec_str.len(), 10);
        assert_eq!(sec_str.hec_code(5), b'C');
    }
}