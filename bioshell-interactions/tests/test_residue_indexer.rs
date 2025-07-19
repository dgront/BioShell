#[cfg(test)]
mod tests {
    use std::io::BufReader;
    use bioshell_interactions::ResidueIndexer;
    use bioshell_pdb::{Deposit, PDBError, Structure};

    const cif_1c5n:  &str = include_str!("./input_files/1c5n.cif");

    #[test]
    fn index_1c5n() -> Result<(), PDBError> {
        let reader = BufReader::new(cif_1c5n.as_bytes());
        let mut strctr = Deposit::from_cif_reader(reader)?.structure()?;
        strctr.remove_ligands();
        let strctr = Structure::from_iterator(&strctr.id_code,
                                              strctr.atoms().iter().filter(|a| a.chain_id == "H").cloned());
        let indexes = ResidueIndexer::from_structure(&strctr);

        Ok(())
    }
}