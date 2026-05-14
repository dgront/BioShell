#[cfg(test)]
mod tests {
    use bioshell_chem::{Atom, BondType, ChemErrors, Molecule};
    use super::*;

    fn atom(name: &str) -> Atom { Atom::neutral(0, 6) }

    fn benzene() -> Result<Molecule, ChemErrors> {
        let mut mol = Molecule::new("benzene");

        for i in 0..6 {
            mol.add_atom(Atom::neutral(i, 6))?;
        }

        for i in 0..6 {
            mol.bind_atoms(i, (i + 1) % 6, BondType::Aromatic)?;
        }

        Ok(mol)
    }

    #[test]
    fn benzene_has_six_atoms_and_six_bonds() -> Result<(), ChemErrors> {
        let mol = benzene()?;

        assert_eq!(mol.molecule_name, "benzene");
        assert_eq!(mol.count_atoms(), 6);
        assert_eq!(mol.count_bonds(), 6);
        Ok(())
    }

    #[test]
    fn benzene_ring_connectivity_is_correct() -> Result<(), ChemErrors> {
        let mol = benzene()?;

        for a in 0..6 {
            let b = (a + 1) % 6;

            assert!(mol.are_bonded(a, b).unwrap());
            assert_eq!(mol.get_bond(a, b).unwrap(), &BondType::Aromatic);
            assert_eq!(mol.count_bonds_for_atom(a), 2);
        }

        assert!(!mol.are_bonded(0, 3)?);
        Ok(())
    }

    #[test]
    fn benzene_adjacency_is_correct() -> Result<(), ChemErrors> {

        let mol = benzene()?;
        let adj = mol.neighbor_indices(0).collect::<Vec<_>>();

        assert_eq!(adj.len(), 2);
        assert!(adj.contains(&1));
        assert!(adj.contains(&5));
        Ok(())
    }

    #[test]
    fn bond_can_be_replaced_and_removed() -> Result<(), ChemErrors> {
        let mut mol = benzene()?;
        mol.add_atom(Atom::neutral(6, 6))?;
        mol.bind_atoms(6, 0, BondType::Single)?;
        assert_eq!(mol.get_bond(6, 0)?, &BondType::Single);

        assert!(mol.remove_bond(6, 0)?);
        assert!(!mol.are_bonded(6, 0)?);
        assert_eq!(mol.count_bonds(), 6);
        Ok(())
    }
}