#[cfg(test)]
mod tests {
    use std::io::BufReader;
    use bioshell_interactions::BackboneHBondMap;
    use bioshell_pdb::Deposit;

    #[allow(non_upper_case_globals)]
    const cif_2fdo:  &str = include_str!("./input_files/2fdo.cif");

    #[test]
    fn hbonds_2fdo() {

        let reader = BufReader::new(cif_2fdo.as_bytes());
        let strctr = Deposit::from_cif_reader(reader).unwrap().structure();
        let hbonds = BackboneHBondMap::new(&strctr);
    }
}