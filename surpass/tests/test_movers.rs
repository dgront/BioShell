use surpass::{HingeMove, MoveProposal, SurpassAlphaSystem};
use bioshell_pdb::calc::dihedral_angle4;
use bioshell_pdb::nerf::restore_linear_chain;

macro_rules! assert_vec3 {
    ($va:expr, $vb:expr, $d:expr) => {
        assert!(($va.x-$vb.x).abs() < $d, "va.x = {}, vb.x = {}", $va.x, $vb.x);
        assert!(($va.y-$vb.y).abs() < $d, "va.y = {}, vb.y = {}", $va.y, $vb.y);
        assert!(($va.z-$vb.z).abs() < $d, "va.z = {}, vb.z = {}", $va.z, $vb.z);
    }
}

macro_rules! assert_delta {
    ($x:expr, $y:expr, $d:expr) => {
        assert!(($x-$y).abs() < $d, "a = {}, b = {}", $x, $y)
    }
}

#[cfg(test)]
mod tests_hinge_move {
    use bioshell_pdb::calc::Vec3;
    use super::*;


    #[test]
    fn make_predefined_hinge_move() {

        make_predefined_hinge_move_N::<1>();
        make_predefined_hinge_move_N::<7>();
        make_predefined_hinge_move_N::<8>();
    }

    fn make_predefined_hinge_move_N<const HINGE_MOVE_SIZE: usize>() {
        // --- a single chain of 10 atoms; HINGE_MOVE_SIZE of them will be moved
        let mut model = SurpassAlphaSystem::new(&[10], 100.0);
        // ---------- Initialize internal coordinates
        let r= [1.0; 10];
        let mut planar = [90.0_f64.to_radians(); 10];
        let mut dihedral = [180.0_f64.to_radians(); 10];
        let mut coords = vec![Vec3::default(); 10];
        restore_linear_chain(&r, &planar, &dihedral, &mut coords);
        for i in 0..10 {
            model.vec3_to_ca(i, &coords[i]);
        }

        let mover: HingeMove<HINGE_MOVE_SIZE> = HingeMove::new(std::f64::consts::PI / 2.0, std::f64::consts::PI / 2.0);
        let mut mp = MoveProposal::new();
        let v0_before = model.ca_to_vec3(0);
        let v1_before = model.ca_to_vec3(1);
        let last_before = model.ca_to_vec3(HINGE_MOVE_SIZE+1);
        mover.compute_move(&model, 1, std::f64::consts::PI / 2.0, &mut mp);
        mp.apply(&mut model);
        let v0_after = model.ca_to_vec3(0);
        let v1_after = model.ca_to_vec3(1);
        let last_after = model.ca_to_vec3(HINGE_MOVE_SIZE+1);

        assert_vec3!(v0_before, v0_after, 0.00001);
        assert_vec3!(last_before, last_after, 0.00001);

        let d = model.ca_to_vec3(1);
        let dihedral_moved = dihedral_angle4(&v1_before, &v0_before, &last_before, &v1_after).to_degrees();
        assert_delta!(dihedral_moved.abs(), 90.0, 0.01);
    }
}