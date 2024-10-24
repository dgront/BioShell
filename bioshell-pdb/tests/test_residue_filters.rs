#[cfg(test)]
mod test_residue_filters {
    use std::io::BufReader;
    use bioshell_cif::read_cif_buffer;
    use bioshell_pdb::residue_filters::{HasAllHeavyAtoms, ResidueFilter};
    use bioshell_pdb::{Deposit, PDBError, ResidueId, Structure};

    #[test]
    fn test_bb_predicate() {
        use bioshell_pdb::residue_filters::{HasCompleteBackbone, ResidueFilter};
        use bioshell_pdb::{PdbAtom, ResidueId, Structure};
        let filter = HasCompleteBackbone;
        let mut strctr = Structure::new("1xyz");
        strctr.push_atom(PdbAtom::from_atom_line("ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N"));
        strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"));
        strctr.push_atom(PdbAtom::from_atom_line("ATOM    516  C   ALA A  69      26.891  29.054  30.649  1.00 15.28           C"));
        strctr.push_atom(PdbAtom::from_atom_line("ATOM    518  CB  ALA A  69      25.155  27.554  29.987  1.00 21.91           C"));
        assert!(!filter.check(&strctr, &ResidueId::new("A", 69, ' ')));
        strctr.push_atom(PdbAtom::from_atom_line("ATOM    517  O   ALA A  69      26.657  29.867  31.341  1.00 20.90           O"));
        assert!(filter.check(&strctr, &ResidueId::new("A", 69, ' ')));
    }

    #[allow(non_upper_case_globals)]
    const cif_2gb1:  &str = include_str!("./test_files/2gb1.cif");


    #[test]
    fn test_heavy_atoms_predicate() -> Result<(), PDBError> {
        let reader = BufReader::new(cif_2gb1.as_bytes());
        let deposit = Deposit::from_cif_reader(reader).unwrap();
        let strctr = deposit.structure();
        let all_heavy = HasAllHeavyAtoms;
        let outcome = all_heavy.check(&strctr, &ResidueId::new("A", 1, ' '));
        assert!(outcome);

        for res_id in strctr.residue_ids() {
            assert!(all_heavy.check(&strctr, res_id));
        }

        Ok(())
    }
}