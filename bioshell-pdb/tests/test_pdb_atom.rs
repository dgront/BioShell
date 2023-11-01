#[cfg(test)]
mod tests {
    use std::string::String;
    use bioshell_pdb::{PdbAtom, SecondaryStructureTypes};

    #[test]
    fn test_new() {
        let atom = PdbAtom::new();
        assert_eq!(atom.serial, 1);
        assert_eq!(atom.name, " CA ");
        assert_eq!(atom.alt_loc, ' ');
        assert_eq!(atom.res_name, "ALA");
        assert_eq!(atom.chain_id, "A");
        assert_eq!(atom.res_seq, 1);
        assert_eq!(atom.i_code, ' ');
        assert_eq!(atom.pos.x, 0.0);
        assert_eq!(atom.pos.y, 0.0);
        assert_eq!(atom.pos.z, 0.0);
        assert_eq!(atom.occupancy, 1.0);
        assert_eq!(atom.temp_factor, 0.0);
        assert_eq!(atom.element, Some(String::from("C")));
        assert_eq!(atom.secondary_struct_type, SecondaryStructureTypes::Coil as u8);
        assert_eq!(atom.charge, None);
    }

    #[test]
    fn test_parse() {
        let atom_line = "ATOM   2831  OE1BGLN A 294C    -27.117  12.343  28.479  1.00  9.58           O  ";
        let atom = PdbAtom::from_atom_line(atom_line);
        assert_eq!(atom.is_hetero_atom,false);
        assert_eq!(atom.serial,2831);
        assert_eq!(atom.element, Some(String::from("O")));
        assert_eq!(atom.alt_loc, 'B');
        assert_eq!(atom.res_name, "GLN");
        assert_eq!(atom.chain_id, "A");
        assert_eq!(atom.res_seq, 294);
        assert_eq!(atom.i_code, 'C');
        assert_eq!(atom.pos.x, -27.117);
        assert_eq!(atom.pos.y, 12.343);
        assert_eq!(atom.pos.z, 28.479);
        assert_eq!(atom.occupancy, 1.0);
        assert_eq!(atom.temp_factor,9.58);
        assert_eq!(atom.charge, None);
        let str = format!("{}",atom);
        assert_eq!(str, atom_line);

        let atom_line = "HETATM 3169  PB  BIS A 499     -10.950  43.372   1.324  1.00 18.00           P  ";
        let atom = PdbAtom::from_atom_line(atom_line);
        let str = format!("{}",atom);
        assert_eq!(str, atom_line);
    }
}