use surpass::{HingeMove, MoveProposal, SurpassAlphaSystem};
use bioshell_pdb::calc::dihedral_angle4;
use bioshell_builder::nerf::restore_linear_chain;

macro_rules! assert_vec3_eq {
    ($va:expr, $vb:expr, $d:expr, $msg:expr) => {
        assert!(($va.x-$vb.x).abs() < $d, "{} : va.x = {}, vb.x = {}", $msg, $va.x, $vb.x);
        assert!(($va.y-$vb.y).abs() < $d, "{} : va.y = {}, vb.y = {}", $msg, $va.y, $vb.y);
        assert!(($va.z-$vb.z).abs() < $d, "{} : va.z = {}, vb.z = {}", $msg, $va.z, $vb.z);
    }
}

macro_rules! assert_vec3_ne {
    ($va:expr, $vb:expr, $d:expr, $msg:expr) => {
        assert!(($va.x-$vb.x).abs() > $d, "{} : va.x = {}, vb.x = {}", $msg, $va.x, $vb.x);
        assert!(($va.y-$vb.y).abs() > $d, "{} : va.y = {}, vb.y = {}", $msg, $va.y, $vb.y);
        assert!(($va.z-$vb.z).abs() > $d, "{} : va.z = {}, vb.z = {}", $msg, $va.z, $vb.z);
    }
}

#[cfg(test)]
mod test_movers {
    use rand::rngs::SmallRng;
    use rand::SeedableRng;
    use bioshell_pdb::calc::Vec3;
    use bioshell_pdb::{assert_delta, load_pdb_file};
    use surpass::{extended_chain, MovedTermini, Mover, TailMove};
    use super::*;


    #[test]
    fn make_predefined_hinge_move() {

        make_predefined_hinge_move_N::<1>();
        make_predefined_hinge_move_N::<7>();
        make_predefined_hinge_move_N::<8>();
    }

    #[allow(non_snake_case)]
    fn make_predefined_hinge_move_N<const HINGE_MOVE_SIZE: usize>() {
        // --- a single chain of 10 atoms; HINGE_MOVE_SIZE of them will be moved
        let mut rnd = SmallRng::seed_from_u64(42);
        let mut model = SurpassAlphaSystem::make_random(&[10], 100.0, &mut rnd);
        // ---------- Initialize internal coordinates
        let r= [1.0; 10];
        let planar = [90.0_f64.to_radians(); 10];
        let dihedral = [180.0_f64.to_radians(); 10];
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

        assert_vec3_eq!(v0_before, v0_after, 0.00001, "first CA");
        assert_vec3_eq!(last_before, last_after, 0.00001, "last CA");

        let dihedral_moved = dihedral_angle4(&v1_before, &v0_before, &last_before, &v1_after).to_degrees();
        assert_delta!(dihedral_moved.abs(), 90.0, 0.01);
    }

    #[test]
    fn make_predefined_tail_move() {
        let mut model = extended_chain(4, 100.0);
        let tail: TailMove<1> = TailMove::new(std::f64::consts::PI / 2.0, std::f64::consts::PI / 2.0);
        let mut tail_prop = MoveProposal::new();
        let n_before = model.ca_to_vec3(0);
        let c_before = model.ca_to_vec3(3);
        tail.compute_move(&mut model, 0, MovedTermini::NTerminal, 1.57, &mut tail_prop);
        tail_prop.apply(&mut model);
        let n_after = model.ca_to_vec3(0);
        let c_after = model.ca_to_vec3(3);
        assert_vec3_ne!(n_before, n_after, 1.5, "N-teminal should move");
        assert_vec3_eq!(c_before, c_after, 0.00001, "C-teminal should NOT move");
        tail.compute_move(&mut model, 0, MovedTermini::CTerminal, 1.57, &mut tail_prop);
        tail_prop.apply(&mut model);
        let c_after = model.ca_to_vec3(3);
        assert_vec3_ne!(c_before, c_after, 1.5, "C-teminal should move");
    }

    #[test]
    fn make_predefined_tail_moves() {
        make_predefined_tail_move_N::<2>();
        make_predefined_tail_move_N::<3>();
    }

    #[allow(non_snake_case)]
    fn make_predefined_tail_move_N<const N_MOVED: usize>() {
        const N_RES: usize = 10;
        let mut model = extended_chain(N_RES, 1000.0);
        let tail: TailMove<N_MOVED> = TailMove::new(std::f64::consts::PI / 2.0, std::f64::consts::PI / 2.0);
        let mut tail_prop = MoveProposal::new();
        let mut backup: MoveProposal<N_MOVED> = MoveProposal::new();
        backup.first_moved_pos = 0;
        backup.backup(&model);
        tail.compute_move(&mut model, 0, MovedTermini::NTerminal, 1.57, &mut tail_prop);
        tail_prop.apply(&mut model);
        for i in 0..N_MOVED {
            assert_ne!(backup.cax[i]+backup.cay[i]+backup.caz[i], model.cax[i]+model.cay[i]+model.caz[i]);
        }
        backup.apply(&mut model);

        backup.first_moved_pos = N_RES-N_MOVED;
        backup.backup(&model);
        tail.compute_move(&mut model, 0, MovedTermini::CTerminal, 1.57, &mut tail_prop);
        tail_prop.apply(&mut model);
        let offset = N_RES-N_MOVED;
        for i in 0..N_MOVED {
            assert_ne!(backup.cax[i]+backup.cay[i]+backup.caz[i], model.cax[i+offset]+model.cay[i+offset]+model.caz[i+offset]);
        }
    }

    #[test]
    fn move_chain() {
        const N: usize = 10;
        let mut rnd = SmallRng::seed_from_u64(42);

        let mut model = SurpassAlphaSystem::make_random(&[N], 1000.0, &mut rnd);
        // ---------- Initialize internal coordinates
        let r= [3.8; N];
        let planar = [120.0_f64.to_radians(); N];
        let dihedral = [180.0_f64.to_radians(); N];
        let mut coords = vec![Vec3::default(); N];
        restore_linear_chain(&r, &planar, &dihedral, &mut coords);
        for i in 0..N {
            model.vec3_to_ca(i, &coords[i]);
        }

        let hinge: HingeMove<4> = HingeMove::new(std::f64::consts::PI / 2.0, std::f64::consts::PI / 2.0);
        let tail: TailMove<1> = TailMove::new(std::f64::consts::PI / 2.0, std::f64::consts::PI / 2.0);
        let mut hinge_prop: MoveProposal<4> = MoveProposal::new();
        let mut tail_prop: MoveProposal<1> = MoveProposal::new();
        // model.to_pdb_file("tra.pdb", false);
        for _i in 0..100000 {
            hinge.propose(&model, &mut rnd, &mut hinge_prop);
            hinge_prop.apply(&mut model);
            tail.propose(&model, &mut rnd, &mut tail_prop);
            tail_prop.apply(&mut model);
            // model.to_pdb_file("tra.pdb", true);
            check_bond_lengths(&model, 3.8);
        }

        for i in 0..N {
            let v = model.ca_to_vec3(i);
            assert_vec3_ne!(v, coords[i], 0.1, format!("atom {} should move",i));
        }
    }

    fn check_bond_lengths(system: &SurpassAlphaSystem, d: f64) {
        for i in 1..system.count_atoms() {
            assert_delta!(system.distance(i-1, i), d, 0.01);
        }
    }
}