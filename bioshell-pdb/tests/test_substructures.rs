#[cfg(test)]
mod test_substructures {
    use std::io::BufReader;
    use bioshell_pdb::{assert_delta, Deposit, PDBError};
    use bioshell_pdb::calc::SubstructureAxis;

    #[allow(non_upper_case_globals)]
    const helix_pdb:  &str = include_str!("./test_files/alphaHelix.pdb");
    #[test]
    fn test_helical_axis() -> Result<(), PDBError> {
        let deposit = Deposit::from_pdb_reader(BufReader::new(helix_pdb.as_bytes()))?;
        let strctr = deposit.structure().unwrap();
        let axis = SubstructureAxis::from_3d_points(&strctr.atoms());
        println!("{} {} {}", axis.begin(), axis.end(), axis.versor());
        Ok(())
    }
}