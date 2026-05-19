#[cfg(test)]
mod tests {
    use bioshell_chem::{BondType, ChemErrors, Element, Molecule};


    #[test]
    fn molecule_from_smiles() -> Result<(), ChemErrors> {
        // let smiles = "c1=ccc=ccc1";
        let smiles = "c1ccccc1";
        let mol = Molecule::from_smiles("test", smiles)?;

        assert_eq!(mol.count_atoms(), 6);
        assert_eq!(mol.count_bonds(), 6);

        for i in 0..6 {
            let atom = mol.get_atom(i).unwrap();
            assert_eq!(atom.element(), Element::C);
            // assert_eq!(atom.charge(), 0);
        }

        for i in 0..6 {
            let j = (i + 1) % 6;
            assert!(mol.are_bonded(i, j).unwrap());
            assert_eq!(mol.get_bond(i, j).unwrap(), &BondType::Aromatic);
        }

        Ok(())
    }

}