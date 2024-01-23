#[cfg(test)]
mod kinematic_tree_tests {
    use std::io::BufReader;
    use bioshell_builder::{InternalCoordinatesDatabase, KinematicAtomTree};
    use bioshell_cif::read_cif_buffer;
    use bioshell_pdb::{assert_delta};
    #[test]
    fn build_backbone() {
        // --- build the database of internal monomers' definitions
        let mut db = InternalCoordinatesDatabase::new();
        let cif_data = read_cif_buffer(BufReader::new(BB_.as_bytes()));
        db.load_from_cif_data(cif_data);
        let cif_data = read_cif_buffer(BufReader::new(CTerm.as_bytes()));
        db.load_from_cif_data(cif_data);

        // --- create a chain builder
        let mut bb_builder = KinematicAtomTree::new();
        let bb_def = db.get_definition("bb_").unwrap();
        bb_builder.add_residue(&bb_def);
        bb_builder.add_residue(&bb_def);
        let cterm_def = db.get_definition("patch_CTerm").unwrap();
        bb_builder.patch_residue(1, &cterm_def);
        let atoms = bb_builder.build_coordinates();
        assert!(atoms.is_ok());
        let atoms = atoms.ok().unwrap();
        assert_eq!(atoms.len(), 9);

        assert_eq!(&(0usize..4), bb_builder.atoms_for_residue(0));
        assert_eq!(&(4usize..9), bb_builder.atoms_for_residue(1));

        // for iatom in 0..9 {
        //     println!("ATOM   {:4} {} GLY {}{:4}    {:8.3}{:8.3}{:8.3}  1.00 99.88           C",
        //              iatom + 1, &bb_builder.atom_name(iatom),
        //              "A", 1, &atoms[iatom].x, &atoms[iatom].y, &atoms[iatom].z);
        // }
    }
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
'bb ' prev ' N  ' prev ' CA ' prev ' C  ' this ' N  ' N  1.328685 114.0  180.0 psi
'bb ' prev ' CA ' prev ' C  ' this ' N  ' this ' CA ' C  1.458001 123.0  180.0 omega
'bb ' prev ' C  ' this ' N  ' this ' CA ' this ' C  ' C  1.523258 110.0 -180.0 phi
'bb ' next ' N  ' this ' CA ' this ' C  ' this ' O  ' O  1.231015 120.0  180.0 -
#
";

    const CTerm: &str = "data_patch_CTerm
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
'CTerm' this ' N  ' this ' CA ' this ' C  ' this ' OXT' O  1.2      116.5  180.0 psi
'CTerm' this ' OXT' this ' CA ' this ' C  ' this ' O  ' O  1.231015 121.0  180.0 -
#
";

}