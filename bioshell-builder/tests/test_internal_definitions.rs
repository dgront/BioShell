#[cfg(test)]
mod internal_atom_definitions_tests {
    use std::io::BufReader;
    use bioshell_builder::{BuilderError, InternalAtomDefinition, InternalCoordinatesDatabase};
    use bioshell_cif::read_cif_buffer;
    use bioshell_io::split_into_strings;
    use bioshell_pdb::assert_delta;

    const BB_: &str = "data_bb_
loop_
_res_name
_atom_a_residue_locator
_atom_a_name
_atom_b_residue_locator
_atom_b_name
_atom_c_residue_locator
_atom_c_name
_atom_d_residue_locator
_atom_d_name
_atom_d_element
_c_d_bond_length
_b_c_d_planar_angle
_a_b_c_d_dihedral_angle
_dihedral_angle_name
'bb ' prev ' N  ' prev ' CA ' prev ' C  ' this ' N  ' N 1.328685 114.0  180.0 psi
'bb ' prev ' CA ' prev ' C  ' prev ' N  ' this ' CA ' C 1.458001 123.0  180.0 omega
'bb ' prev ' C  ' this ' N  ' this ' CA ' this ' C  ' C 1.523258 110.0 -180.0 phi
'bb ' this ' N  ' this ' CA ' this ' C  ' next ' N  ' N 1.328685 114.0  180.0 psi
'bb ' next ' N  ' this ' CA ' this ' C  ' this ' O  ' O 1.231015 121.0  180.0 -
#
";

    #[test]
    fn test_cif_loading() -> Result<(), BuilderError> {
        let mut db = InternalCoordinatesDatabase::new();
        let cif_data = read_cif_buffer(BufReader::new(BB_.as_bytes()))?;
        db.load_from_cif_data(cif_data)?;
        let bb_def = db.get_definition("bb_");
        assert!(bb_def.is_some());
        let bb_def = bb_def.unwrap();
        assert_eq!(bb_def.len(), 5);
        Ok(())
    }

    #[test]
    fn test_internal_definitions() -> Result<(), BuilderError> {
        let data: Vec<String> = split_into_strings("'ALA' this ' N  ' this ' CA ' this ' C  ' next ' N  ' N 1.328685 114.0  180.0 psi", true);
        let data_str: Vec<&str> = data.iter().map(AsRef::as_ref).collect();
        let def = InternalAtomDefinition::from_strings(&data_str)?;
        assert_eq!(def.a_name, " N  ".to_string());
        assert_eq!(def.name, " N  ".to_string());
        assert_delta!(def.r, 1.328685, 0.000001);
        assert_delta!(def.planar, 114.0_f64.to_radians(), 0.0001);
        assert_delta!(def.dihedral, 180.0_f64.to_radians(), 0.0001);

        let data: Vec<String> = split_into_strings("'bb ' next ' N  ' this ' CA ' this ' C  ' this ' O  ' O 1.231015 121.0  180.0 -", true);
        let data_str: Vec<&str> = data.iter().map(AsRef::as_ref).collect();
        let def = InternalAtomDefinition::from_strings(&data_str)?;
        assert_eq!(def.a_name, " N  ".to_string());
        assert_eq!(def.name, " O  ".to_string());
        assert_delta!(def.r, 1.231015, 0.00001);
        assert_delta!(def.planar, 121.0_f64.to_radians(), 0.0001);
        assert_delta!(def.dihedral, 180.0_f64.to_radians(), 0.0001);
        Ok(())
    }
}