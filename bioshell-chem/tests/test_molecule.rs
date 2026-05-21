#[cfg(test)]
mod tests {
    use bioshell_chem::{Atom, BondType, ChemErrors, Element, load_molecule, Molecule, write_mol2};
    use bioshell_chem::icoords::{ZMatrix};
    use bioshell_core::assert_delta;
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
        mol.add_atom(Atom::neutral(6, Element::C))?;
        mol.bind_atoms(6, 0, BondType::Single)?;
        assert_eq!(mol.get_bond(6, 0)?, &BondType::Single);

        assert!(mol.remove_bond(6, 0)?);
        assert!(!mol.are_bonded(6, 0)?);
        assert_eq!(mol.count_bonds(), 6);
        Ok(())
    }

    #[test]
    fn molecule_from_cif() -> Result<(), ChemErrors> {

        let reader = open_file("tests/test_files/MBN.cif")?;
        let mol = bioshell_chem::molecule_from_cif(reader)?;
        assert_eq!(mol.molecule_name, "TOLUENE");
        assert_eq!(mol.code, Some("MBN".to_string()));
        assert_eq!(mol.count_atoms(), 15);
        assert_eq!(mol.count_bonds(), 15);
        Ok(())
    }

    #[test]
    fn molecule_from_sdf() -> Result<(), ChemErrors> {

        let reader = open_file("tests/test_files/toluene.sdf")?;
        let mol = bioshell_chem::molecule_from_sdf(reader)?;
        assert_eq!(mol.count_atoms(), 15);
        assert_eq!(mol.count_bonds(), 15);
        Ok(())
    }

    /// The TOL.gro file content
    const TOL_GRO: &str = r#"LIGPARGEN GENERATED GRO FILE
   15
    1TOL    C00    1  -0.208   0.000  -0.000
    1TOL    H01    2  -0.316   0.000  -0.000
    1TOL    C02    3  -0.139   0.121  -0.000
    1TOL    H03    4  -0.193   0.215  -0.000
    1TOL    C04    5  -0.139  -0.121  -0.000
    1TOL    H05    6  -0.193  -0.215  -0.000
    1TOL    C06    7   0.002   0.121  -0.000
    1TOL    H07    8   0.055   0.215   0.000
    1TOL    C08    9   0.002  -0.121  -0.000
    1TOL    H09   10   0.055  -0.215   0.000
    1TOL    C0A   11   0.072  -0.000  -0.000
    1TOL    C0B   12   0.222  -0.000  -0.001
    1TOL    H0C   13   0.258   0.001   0.105
    1TOL    H0D   14   0.260   0.090  -0.054
    1TOL    H0E   15   0.260  -0.091  -0.052
   1.00000   1.00000   1.00000
"#;
    #[test]
    fn molecule_from_itp() -> Result<(), ChemErrors> {

        let reader = open_file("tests/test_files/TOL.itp")?;
        let mut mol = bioshell_chem::molecule_from_itp(reader)?;
        assert_eq!(mol.count_atoms(), 15);
        assert_eq!(mol.count_bonds(), 15);

        let lines: Vec<&str> =  TOL_GRO.lines().collect();
        for (i,line) in lines[2..17].iter().enumerate() {
            let f: Vec<&str> = line.split_whitespace().collect();
            let x = f[3].parse::<f64>().unwrap() * 10.0; // convert from nm to A
            let y = f[4].parse::<f64>().unwrap() * 10.0;
            let z = f[5].parse::<f64>().unwrap() * 10.0;
            mol.set_pos3(i, x, y, z);
        }
        assert_delta!(mol.pos(13).x, 2.60, 0.001);
        assert_delta!(mol.pos(13).y, 0.90, 0.001);
        assert_delta!(mol.pos(13).z, -0.54, 0.001);
        let _zmat = ZMatrix::from_molecule(&mut mol, 0, 2, 4)?;
        // println!("{}", &_zmat);

        Ok(())
    }

    #[test]
    fn molecule_to_mol2() -> Result<(), ChemErrors> {
        let expected = "@<TRIPOS>MOLECULE
ETHANOL
    9     8 0 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
      1 C1        15.2120    49.1980     7.4910 C         1 MOL
      2 C2        16.0690    50.3860     7.1040 C         1 MOL
      3 O3        15.8610    48.1850     8.2560 O         1 MOL
      4 H4        14.3750    49.5790     8.0940 H         1 MOL
      5 H5        14.8580    48.7310     6.5600 H         1 MOL
      6 H6        15.4670    51.0980     6.5200 H         1 MOL
      7 H7        16.4420    50.8800     8.0130 H         1 MOL
      8 H8        16.9200    50.0420     6.4980 H         1 MOL
      9 H9        15.2440    47.4880     8.4470 H         1 MOL
@<TRIPOS>BOND
     1     1     2 1
     2     1     3 1
     3     1     4 1
     4     1     5 1
     5     2     6 1
     6     2     7 1
     7     2     8 1
     8     3     9 1";

        let mut buf = Vec::new();
        let et_oh = load_molecule("./tests/test_files/EOH.cif")?;
        write_mol2(&et_oh, &mut buf)?;
        let text = String::from_utf8(buf).unwrap();

        for (line_no, (actual_line, expected_line)) in text.lines().zip(expected.lines()).enumerate() {
            assert_eq!(
                actual_line.trim(),
                expected_line.trim(),
                "Mismatch at line {}",
                line_no + 1
            );
        }

        Ok(())
    }
}