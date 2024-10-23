#[cfg(test)]
mod tests {
    use std::io::BufReader;
    use itertools::assert_equal;
    use bioshell_pdb::{load_cif_reader, load_pdb_reader, PdbAtom, Structure};
    use bioshell_pdb::pdb_atom_filters::{ByChain, ByResidue, ByResidueType, IsAromatic, IsBackbone, IsCA, MatchAll, MatchAny, PdbAtomPredicate, PdbAtomPredicate2, SameResidue};

    #[allow(non_upper_case_globals)]
    const pdb_2gb1:  &str = include_str!("./test_files/2gb1.pdb");
    #[allow(non_upper_case_globals)]
    const pdb_2fdo:  &str = include_str!("./test_files/2fdo.cif");

    #[allow(non_upper_case_globals)]
    const lines: [&str; 15] = [
        "ATOM    514  N   ALA A  68      26.532  28.200  28.365  1.00 17.85           N",
        "ATOM    515  CA  ALA A  68      25.790  28.757  29.513  1.00 16.12           C",
        "ATOM    516  C   ALA A  68      26.891  29.054  30.649  1.00 15.28           C",
        "ATOM    517  O   ALA A  68      26.657  29.867  31.341  1.00 20.90           O",
        "ATOM    518  CB  ALA A  68      25.155  27.554  29.987  1.00 21.91           C",
        "ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N",
        "ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C",
        "ATOM    516  C   ALA A  69      26.891  29.054  30.649  1.00 15.28           C",
        "ATOM    517  O   ALA A  69      26.657  29.867  31.341  1.00 20.90           O",
        "ATOM    518  CB  ALA A  69      25.155  27.554  29.987  1.00 21.91           C",
        "ATOM    514  N   ALA B  68      26.532  28.200  28.365  1.00 17.85           N",
        "ATOM    515  CA  ALA B  68      25.790  28.757  29.513  1.00 16.12           C",
        "ATOM    516  C   ALA B  68      26.891  29.054  30.649  1.00 15.28           C",
        "ATOM    517  O   ALA B  68      26.657  29.867  31.341  1.00 20.90           O",
        "ATOM    518  CB  ALA B  68      25.155  27.554  29.987  1.00 21.91           C"];

    #[test]
    fn test_backbone_filtering() {
        let atoms: Vec<PdbAtom> = lines.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
        let mut strctr = Structure::from_iterator("1xyz", atoms.iter());
        let bb = IsBackbone {};
        let bb_count = strctr.atoms().iter().filter(|a| bb.check(a)).count();
        assert_eq!(bb_count, 12);
    }

    #[test]
    fn structure_from_iterator() {
        let atoms: Vec<PdbAtom> = lines.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
        let strctr = Structure::from_iterator("1xyz", atoms.iter());
        let bb = IsBackbone {};
        let bb_strctr = Structure::from_iterator("1xyz", strctr.atoms().iter().filter(|a| bb.check(a)));
        assert_eq!(bb_strctr.count_atoms(), 12);
    }


    #[test]
    fn iterate_over_residues() {
        let atoms: Vec<PdbAtom> = lines.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
        let strctr = Structure::from_iterator("1xyz", atoms.iter());

        let same_res = SameResidue {};
        let n_res = strctr.atoms().windows(2).filter(|a| !same_res.check(&a[0], &a[1])).count() + 1;
        assert_eq!(n_res, 3);
    }

    #[test]
    fn residues_of_a_chain() {
        let atoms: Vec<PdbAtom> = lines.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
        let strctr = Structure::from_iterator("1xyz", atoms.iter());

        let chain_a = ByChain::new("A");
        let residues = Structure::residue_ids_from_atoms(strctr.atoms().iter().filter(|a| chain_a.check(a)));
        assert_eq!(residues.len(), 2);

        for res_idx in residues {
            let selector = ByResidue::new(res_idx);
            assert_eq!(strctr.atoms().iter().filter(|&a| selector.check(a)).count(), 5);
        }

        let chain_b = ByChain::new("B");
        let residues = Structure::residue_ids_from_atoms(strctr.atoms().iter().filter(|a| chain_b.check(a)));
        assert_eq!(residues.len(), 1);

        for res_idx in residues {
            let selector = ByResidue::new(res_idx);
            assert_eq!(strctr.atoms().iter().filter(|&a| selector.check(a)).count(), 5);
        }
    }

    #[test]
    fn test_chain_filtering() {
        let by_chain_a = ByChain::new("A");
        let strctr = load_cif_reader(BufReader::new(pdb_2fdo.as_bytes())).unwrap();
        let chain_a_atoms = strctr.atoms().iter().filter(|a| by_chain_a.check(a)).collect::<Vec<_>>();
        assert_eq!(chain_a_atoms.len(), 730);
    }

    #[test]
    fn aromatic_residues() {
        let strctr = load_pdb_reader(BufReader::new(pdb_2gb1.as_bytes())).unwrap();

        let mut filter = MatchAll::new();
        filter.add_predicate(Box::new(IsAromatic));
        filter.add_predicate(Box::new(IsCA));

        let aro_res_count: usize = strctr.atoms().iter().filter(|a| filter.check(a)).count();
        assert_eq!(aro_res_count, 6);
    }
    #[test]
    fn ala_or_gly() {
        let strctr = load_pdb_reader(BufReader::new(pdb_2gb1.as_bytes())).unwrap();

        let mut filter = MatchAny::new();
        filter.add_predicate(Box::new(ByResidueType::new("ALA")));
        filter.add_predicate(Box::new(ByResidueType::new("GLY")));

        let ala_gly_res_count: usize = strctr.atoms().iter().filter(|a| filter.check(a)).count();
        assert_eq!(ala_gly_res_count, 88);
    }
}