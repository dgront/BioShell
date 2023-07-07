#[cfg(test)]
mod tests {
    use bioshell_numerical::Vec3;
    use bioshell_pdb::pdb_line_parser::PdbLineParser;
    use bioshell_pdb::PdbAtom;

    #[test]
    fn test_parse_atom() {
        let pdb_line = "ATOM   2831  OE1 GLN A 294     -27.117  12.343  28.479  1.00  9.58           O  ";
        let elements = PdbLineParser::parse_atom(pdb_line).unwrap();
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
        let elements = PdbLineParser::parse_atom(pdb_line);
        assert_eq!(elements, None);
    }

    #[test]
    fn test_assemble_atom() {
        let mut pdb_atom = PdbAtom::new();
        pdb_atom.atom_serial_no = Some(2831);
        pdb_atom.atom_symbol = "OE1".to_string();
        pdb_atom.residue_name = "GLN".to_string();
        pdb_atom.chain_name = "A".to_string();
        pdb_atom.residue_no = Some(294);
        pdb_atom.coordinate = Vec3::new(-27.117, 12.343, 28.479);
        pdb_atom.occupancy = Some(1.0);
        pdb_atom.temperature_factor = Some(9.58);
        pdb_atom.charge_of_the_atom = "O".to_string();

        let pdb_line = PdbLineParser::assemble_atom(&pdb_atom);
        assert_eq!(pdb_line, "ATOM     2831 OE1  GLN A 294      -27.117  12.343  28.479  1.00  9.58          O ");
    }
}