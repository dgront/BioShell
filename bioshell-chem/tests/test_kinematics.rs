#[cfg(test)]
mod tests {
    use bioshell_chem::{Atom, BondType, ChemErrors, Element, load_molecule, Molecule, write_mol2};
    use bioshell_chem::icoords::{KinematicAtomTree, ZMatrix};
    use bioshell_core::io::open_file;
    use bioshell_core::{assert_delta, assert_vec3_eq, Vec3};

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

    #[test]
    fn reconstruct_cartesians() -> Result<(), ChemErrors> {

        // --- read in a molecule
        // Atom order in the CIF file is: C1, C2, O and then hydrogens
        // where the chain goes by: C2-C1-O-H i.e. O is connected to C1
        let mut EtOH = load_molecule("./tests/test_files/EOH.cif")?;
        // --- distance between methyl C2 and hydroxyl H
        let ch_dist = EtOH.pos(1).distance_to(EtOH.pos(8));
        assert_delta!(ch_dist, 3.299, 0.001);

        // --- clone its Cartesians
        let xyz_clone = EtOH.positions().to_vec();

        // --- compute internals
        let mut zmatrix = ZMatrix::from_molecule(&mut EtOH, 2, 0, 1)?;
        let inc_clone = zmatrix.internal_coordinates().to_vec();

        // --- zero the molecule's Cartesians
        let mut zeros = vec![Vec3::default();xyz_clone.len()];
        zeros[0].set(&xyz_clone[0]);
        for (i, p) in zeros.iter().enumerate() {
            EtOH.set_pos3(i,p.x, p.y, p.z);
        }
        // --- test if they are really zero
        let zero = Vec3::default();
        for ai in &EtOH.positions()[1..] {
            assert_vec3_eq!(ai, zero, 0.00001, "atom is not zero");
        }

        // --- create a new z-matrix with the zero-ed molecule
        let mut zmatrix = ZMatrix::from_molecule(&mut EtOH, 2, 0, 1)?;
        // --- set the old internals
        zmatrix.set_internals(&inc_clone)?;

        let ch_dist_again = EtOH.pos(1).distance_to(EtOH.pos(8));
        assert_delta!(ch_dist, ch_dist_again, 0.0001);
        let c = EtOH.get_atom(0)?;
        let h = EtOH.get_atom(8)?;
        let mut buf = Vec::new();
        write_mol2(&EtOH, &mut buf);
        let text = String::from_utf8(buf).unwrap();
        println!("{}", text);

        Ok(())
    }
}