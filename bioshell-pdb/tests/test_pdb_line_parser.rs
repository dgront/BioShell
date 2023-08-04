#[cfg(test)]
mod tests {
    use bioshell_numerical::Vec3;
    use bioshell_pdb::pdb_atom_line_parser::AtomLineParser;
    use bioshell_pdb::PdbAtom;

    #[test]
    fn test_parse_hetero_atom() {
        let pdb_line = "HETATM 3169  PB  BIS A 499     -10.950  43.372   1.324  1.00 18.00           P  ";
        let elements = AtomLineParser::parse_atom(pdb_line).unwrap();
        assert_eq!(elements[0], "HETATM");
        assert_eq!(elements[1], "3169");
        assert_eq!(elements[2], "PB");
        assert_eq!(elements[3], "");
        assert_eq!(elements[4], "BIS");
        assert_eq!(elements[5], "A");
        assert_eq!(elements[6], "499");
        assert_eq!(elements[7], "");
        assert_eq!(elements[8], "-10.950");
        assert_eq!(elements[9], "43.372");
        assert_eq!(elements[10], "1.324");
        assert_eq!(elements[11], "1.00");
        assert_eq!(elements[12], "18.00");
        assert_eq!(elements[13], "");
        assert_eq!(elements[14], "P");
        assert_eq!(elements[15], "");


    }

    #[test]
    fn test_parse_atom() {
        let pdb_line = "ATOM   2831  OE1 GLN A 294     -27.117  12.343  28.479  1.00  9.58           O  ";
        let elements = AtomLineParser::parse_atom(pdb_line).unwrap();
        assert_eq!(elements[0], "ATOM");
        assert_eq!(elements[1], "2831");
        assert_eq!(elements[2], "OE1");
        assert_eq!(elements[3], "");
        assert_eq!(elements[4], "GLN");
        assert_eq!(elements[5], "A");
        assert_eq!(elements[6], "294");
        assert_eq!(elements[7], "");
        assert_eq!(elements[8], "-27.117");
        assert_eq!(elements[9], "12.343");
        assert_eq!(elements[10], "28.479");
        assert_eq!(elements[11], "1.00");
        assert_eq!(elements[12], "9.58");
        assert_eq!(elements[13], "");
        assert_eq!(elements[14], "O");
        assert_eq!(elements[15], "");

        // Test with a line that doesn't start with "ATOM"
        let pdb_line = "HETATM 2831  OE1 GLN A 294     -27.117  12.343  28.479  1.00  9.58           O  ";
        let elements = AtomLineParser::parse_atom(pdb_line).unwrap();
        assert_eq!(elements, ["HETATM", "2831", "OE1", "", "GLN", "A", "294", "", "-27.117", "12.343", "28.479", "1.00", "9.58", "", "O", ""]);
    }

    #[test]
    fn test_assemble_atom() {
        let mut pdb_atom = PdbAtom::new();
        pdb_atom.atom_serial_no_1 = Some(2831);
        pdb_atom.atom_symbol_2 = "OE1".to_string();
        pdb_atom.residue_name_7 = "GLN".to_string();
        pdb_atom.chain_name_8 = "A".to_string();
        pdb_atom.residue_no_9 = Some(294);
        pdb_atom.coordinate_11 = Vec3::new(-27.117, 12.343, 28.479);
        pdb_atom.occupancy_12 = Some(1.0);
        pdb_atom.temperature_factor_13 = Some(9.58);
        pdb_atom.charge_of_the_atom_16 = "O".to_string();

        let pdb_line = AtomLineParser::assemble_atom(&pdb_atom);
        assert_eq!(pdb_line, "ATOM   2831 OE1  GLN A 294  -27.117  12.343  28.479  1.00  9.58            O ");
    }
}