#[cfg(test)]
mod tests {
    use bioshell_chem::{Atom, BondType, ChemErrors, Element, load_molecule, Molecule};
    use bioshell_chem::icoords::{KinematicAtomTree, ZMatrix};
    use bioshell_core::io::open_file;

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
        let benzene_tree = KinematicAtomTree::from_molecule(&mol, 0, 1, 2)?;
        for idef in benzene_tree.atoms() {
            let atom = mol.get_atom(idef.atom)?;
            eprintln!("{:3} {}: {:3} {:3} {:3}", idef.atom, atom.element(),idef.a, idef.b, idef.c);
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

    #[test]
    fn create_z_matrix() -> Result<(), ChemErrors> {
        let mut mol = load_molecule("tests/test_files/MBN.cif")?;
        let _zmat = ZMatrix::from_molecule(&mut mol, 0, 1, 2)?;
        Ok(())
    }
}