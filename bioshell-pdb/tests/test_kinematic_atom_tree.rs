#[cfg(test)]
mod nerf_test {
    use std::f64::consts::PI;
    use bioshell_pdb::calc::{Vec3, dihedral_angle4, planar_angle3};
    // use bioshell_pdb::nerf::{Kine};
    use bioshell_pdb::{assert_delta, PdbAtom};
    use bioshell_pdb::nerf::{InternalAtomDefinition, KinematicAtomTree};

    #[test]
    fn build_backbone() {
        let mut bb_builder = KinematicAtomTree::new();
        let phi = -57.8_f64.to_radians(); // -180.0_f64.to_radians();
        let psi = -47.0_f64.to_radians(); // -180.0_f64.to_radians();
        let omega = 180.0_f64.to_radians();
        let n_def = InternalAtomDefinition::from_properties(
            " N  ", " N  ", " CA ", " C  ",
            1.328685, 114.0_f64.to_radians(), psi);
        let ca_def = InternalAtomDefinition::from_properties(
            " CA ", " CA ", " C  ", " N  ",
            1.458001, 123.0_f64.to_radians(), omega);
        let c_def = InternalAtomDefinition::from_properties(
            " C  ", " C  ", " N  ", " CA ",
            1.523258, 110.0_f64.to_radians(), phi);
        let o_def = InternalAtomDefinition::from_properties(
            " O  ", " N  ", " CA ", " C  ",
            1.231015, 121.0_f64.to_radians(), 180.0_f64.to_radians());

        bb_builder.add_atom(&n_def,0);
        let mut _elements = vec!["N "];
        let _elements4 = ["C ", "C ", "O ", "N "];
        for i in 0..10 {
            bb_builder.add_atom(&ca_def,i*4 + 1);
            bb_builder.add_atom(&c_def,i*4 + 2);
            bb_builder.add_atom(&n_def,i*4 + 4);
            bb_builder.add_atom(&o_def,i*4 + 3);
            _elements.extend_from_slice(&_elements4);
        }

        let chain = bb_builder.restore_atoms();
        // for iatom in 0..chain.len() {
        //     let i_residue = iatom / 4 + 1;
        //     println!("ATOM   {:4} {} GLY {}{:4}    {:8.3}{:8.3}{:8.3}  1.00 99.88           {:}",
        //             iatom + 1, bb_builder.name(iatom).unwrap(), "A", i_residue, &chain[iatom].x, &chain[iatom].y, &chain[iatom].z, &_elements[iatom]);
        // }
        let n_n_distance = chain[0].distance_to(&chain.last().unwrap());
        assert_delta!(n_n_distance, 14.99, 0.001);

        let phi = -180.0_f64.to_radians();
        let psi = 180.0_f64.to_radians();
        for i in 0..10 {
            bb_builder.set_dihedral(i*4 + 3, psi);
            bb_builder.set_dihedral(i*4 + 2, phi);
        }
        let chain = bb_builder.restore_atoms();
        // for iatom in 0..chain.len() {
        //     let i_residue = iatom / 4 + 1;
        //     println!("ATOM   {:4} {} GLY {}{:4}    {:8.3}{:8.3}{:8.3}  1.00 99.88           {:}",
        //             iatom + 1, bb_builder.name(iatom).unwrap(), "A", i_residue, &chain[iatom].x, &chain[iatom].y, &chain[iatom].z, &_elements[iatom]);
        // }
        let n_n_distance = chain[0].distance_to(&chain.last().unwrap());
        assert_delta!(n_n_distance, 36.207, 0.001);
    }
}