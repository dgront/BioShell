#[cfg(test)]
mod tests {
    use bioshell_numerical::Vec3;
    use bioshell_pdb::PdbAtom;
    use super::*;

    #[test]
    fn test_new() {
        let atom = PdbAtom::new();
        assert_eq!(atom.atom_serial_no, Some(0));
        assert_eq!(atom.atom_symbol, "");
        assert_eq!(atom.atom_position, "");
        assert_eq!(atom.atom_no_in_the_branch, Some(0));
        assert_eq!(atom.connected_to_atom_no_in_the_branch, Some(0));
        assert_eq!(atom.alt_loc_indicator, "");
        assert_eq!(atom.residue_name, "");
        assert_eq!(atom.chain_name, "");
        assert_eq!(atom.residue_no, Some(0));
        assert_eq!(atom.insertion_code, "");
        assert_eq!(atom.coordinate.x, 0.0);
        assert_eq!(atom.coordinate.y, 0.0);
        assert_eq!(atom.coordinate.z, 0.0);
        assert_eq!(atom.occupancy, Some(0.0));
        assert_eq!(atom.temperature_factor, Some(0.0));
        assert_eq!(atom.segment_identifier, "");
        assert_eq!(atom.segment_identifier_symbol, "");
        assert_eq!(atom.charge_of_the_atom, "");
    }

    #[test]
    fn test_parse() {
        let atom_line = "ATOM   2831  OE1 GLN A 294     -27.117  12.343  28.479  1.00  9.58           O  ";
        let atom = PdbAtom::parse(atom_line);
        assert_eq!(atom.atom_serial_no, Some(2831));
        assert_eq!(atom.atom_symbol, "O");
        assert_eq!(atom.atom_position, "E");
        assert_eq!(atom.atom_no_in_the_branch, Some(1));
        assert_eq!(atom.connected_to_atom_no_in_the_branch, None);
        assert_eq!(atom.alt_loc_indicator, "");
        assert_eq!(atom.residue_name, "GLN");
        assert_eq!(atom.chain_name, "A");
        assert_eq!(atom.residue_no, Some(294));
        assert_eq!(atom.insertion_code, "");
        assert_eq!(atom.coordinate.x, -27.117);
        assert_eq!(atom.coordinate.y, 12.343);
        assert_eq!(atom.coordinate.z, 28.479);
        assert_eq!(atom.occupancy, Some(1.0));
        assert_eq!(atom.temperature_factor, Some(9.58));
        assert_eq!(atom.segment_identifier, "");
        assert_eq!(atom.segment_identifier_symbol, "");
        assert_eq!(atom.charge_of_the_atom, "O");
    }

    #[test]
    fn test_to_csv_string() {
        let mut atom = PdbAtom::new();
        atom.atom_serial_no = Some(2831);
        atom.atom_symbol = "O".to_string();
        atom.atom_position = "E".to_string();
        atom.atom_no_in_the_branch = Some(1);
        atom.residue_name = "GLN".to_string();
        atom.chain_name = "A".to_string();
        atom.residue_no = Some(294);
        atom.coordinate = Vec3::new(-27.117, 12.343, 28.479);
        atom.occupancy = Some(1.0);
        atom.temperature_factor = Some(9.58);
        atom.charge_of_the_atom = "O".to_string();
        let csv_string = atom.to_csv_string();
        assert_eq!(csv_string, "2831,O,E,1,0,,GLN,A,294,,,-27.117,12.343,28.479,1,9.58,,O");
    }

    #[test]
    fn test_header() {
        let header = PdbAtom::header();
        assert_eq!(header, "AtomSerialNo,AtomSymbol,AtomPosition,AtomNoInTheBranch,ConnectedToAtomNoInTheBranch,AltLocIndicator,ResidueName,ChainName,ResidueNo,InsertionCode,Coordinate.X,Coordinate.Y,Coordinate.Z,Occupancy,TemperatureFactor,SegmentIdentifier,ChargeOfTheAtom");
    }
}