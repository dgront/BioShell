#[cfg(test)]
mod tests {
    use bioshell_chem::{Atom, BondType, ChemErrors, Element, Molecule};
    use bioshell_chem::icoords::KinematicAtomChain;
    use bioshell_io::open_file;

    fn benzene() -> Result<Molecule, ChemErrors> {
        let mut mol = Molecule::new("benzene");

        for i in 0..6 {
            mol.add_atom(Atom::neutral(i, Element::C))?;
        }

        for i in 0..6 {
            mol.bind_atoms(i, (i + 1) % 6, BondType::Aromatic)?;
        }

        Ok(mol)
    }

    #[test]
    fn benzene_kinematics() -> Result<(), ChemErrors> {
        let mut mol = benzene()?;
        mol.add_hydrogens()?;
        let benzene_tree = KinematicAtomChain::from_molecule(&mol, 0, 1, 2)?;
        for idef in benzene_tree.atoms() {
            let atom = mol.get_atom(idef.atom).ok_or_else(|| ChemErrors::InvalidAtomIndex(idef.atom))?;
            eprintln!("{:3} {}: {:3} {:3} {:3}", idef.atom, atom.element(),idef.i, idef.j, idef.k);
        }
        Ok(())
    }


    #[test]
    fn molecule_from_cif() -> Result<(), ChemErrors> {

        let reader = open_file("tests/test_files/CLR.cif")?;
        let mol = bioshell_chem::molecule_from_cif(reader)?;
        assert_eq!(mol.molecule_name, "CHOLESTEROL");
        assert_eq!(mol.code, Some("CLR".to_string()));
        assert_eq!(mol.count_atoms(), 74);
        assert_eq!(mol.count_bonds(), 77);
        Ok(())
    }
}