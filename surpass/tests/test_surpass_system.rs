use std::io::BufReader;
use std::string::String;
use rand::rngs::SmallRng;
use rand::SeedableRng;
use bioshell_pdb::calc::{planar_angle3, Vec3};
use bioshell_pdb::load_pdb_reader;
use bioshell_pdb::nerf::restore_linear_chain;
use surpass::{extended_chain, SurpassAlphaSystem};

#[allow(non_upper_case_globals)]
const pdb_txt: &str = "ATOM      2  CA  MET A   1     -13.296   0.028   3.924  1.00  0.43           C
ATOM     21  CA  THR A   2      -9.669  -0.447   4.998  1.00  0.19           C
ATOM     35  CA  TYR A   3      -7.173  -2.314   2.811  1.00  0.08           C
ATOM     56  CA  LYS A   4      -3.922  -3.881   4.044  1.00  0.10           C
ATOM     78  CA  LEU A   5      -0.651  -2.752   2.466  1.00  0.11           C
ATOM     97  CA  ILE A   6       2.338  -5.105   2.255  1.00  0.13           C
ATOM      2  CA  MET B   1     -13.296   0.028   3.924  1.00  0.43           C
ATOM     21  CA  THR B   2      -9.669  -0.447   4.998  1.00  0.19           C
ATOM     35  CA  TYR B   3      -7.173  -2.314   2.811  1.00  0.08           C
ATOM     56  CA  LYS B   4      -3.922  -3.881   4.044  1.00  0.10           C
ATOM     78  CA  LEU B   5      -0.651  -2.752   2.466  1.00  0.11           C";

macro_rules! assert_delta {
    ($x:expr, $y:expr, $d:expr) => {
        assert!(($x-$y).abs() < $d, "a = {}, b = {}", $x, $y)
    }
}

#[test]
fn system_from_pdb() {

    let strctr = load_pdb_reader(BufReader::new(pdb_txt.as_bytes())).unwrap();

    let model = SurpassAlphaSystem::from_pdb_structure(&strctr, 100.0);
    // --- check if chains are properly assigned
    assert_eq!(model.count_chains(), 2);
    assert_eq!(model.chain(7), 1);
    assert_eq!(*model.chain_id(1).unwrap(), String::from("B"));

    // --- check coordinates
    assert_delta!(model.int_to_real(model.bbx[0]), -13.296, 0.001);
    assert_delta!(model.int_to_real(model.bby[0]), 0.028, 0.001);
    assert_delta!(model.int_to_real(model.bbz[0]), 3.924, 0.001);
    assert_delta!(model.int_to_real(model.bbx[10]), -0.651, 0.001);
    assert_delta!(model.int_to_real(model.bby[10]), -2.752, 0.001);
    assert_delta!(model.int_to_real(model.bbz[10]), 2.466, 0.001);

    // --- compute the distance between two atoms
    assert_delta!(model.distance(0, 1), 3.812, 0.01);
}

/// Make a system of 4 chains and check whether atoms are correctly assigned to chains
#[test]
fn test_4_chains() {

    let mut rnd = SmallRng::from_entropy();
    let model = SurpassAlphaSystem::make_random(&[10, 10, 10, 10], 100.0, &mut rnd);
    for i in 0..4 {
        assert_eq!(model.chain_residues(i).start, i*10);
        assert_eq!(model.chain_residues(i).end, (i+1)*10);
        for j in i*10..(i+1)*10 {
            assert_eq!(model.chain(j), i as u16);
        }
    }
    // model.to_pdb_file("s.pdb", false);
}

#[test]
fn test_extended_chain() {
    let mut rnd = SmallRng::from_entropy();
    let model = extended_chain(10, 100.0);
    assert_eq!(model.chain_residues(0).start, 0);
    assert_eq!(model.chain_residues(0).end, 10);
    for j in 0..10 {
        assert_eq!(model.chain(j), 0);
    }
}

#[test]
fn test_coords_operations() {
    let mut model = SurpassAlphaSystem::new(&[3], 100.0);

    model.bbx[2] = model.real_to_int(6.0);
    assert_delta!(model.int_to_real(model.bbx[2]), 6.0, 0.00001);
    // --- set X coordinate beyond the box
    model.bbx[2] = model.real_to_int(51.0);
    assert_delta!(model.int_to_real(model.bbx[2]), -49.0, 0.00001);
    // --- set X coordinate inside the box but close to the negative end
    model.bbx[2] = model.real_to_int(-49.0);
    assert_delta!(model.int_to_real(model.bbx[2]), -49.0, 0.00001);
    // --- set X coordinate outside the box on the negative side
    model.bbx[2] = model.real_to_int(-51.0);
    assert_delta!(model.int_to_real(model.bbx[2]), 49.0, 0.00001);
}

#[test]
fn test_diatance_evaluation() {
    let mut model = SurpassAlphaSystem::new(&[3], 100.0);
    for i in 0..3 {
        model.bbx[i] = model.real_to_int(0.0);
        model.bby[i] = model.real_to_int(0.0);
        model.bbz[i] = model.real_to_int(0.0);
    }
    model.bbx[2] = model.real_to_int(6.0);
    model.bbx[1] = model.real_to_int(3.0);
    model.bby[1] = model.real_to_int(4.0);
    assert_delta!(model.distance(0,1), 5.0, 0.00001);
    assert_delta!(model.distance(2,1), 5.0, 0.00001);
}

#[test]
fn test_bond_adjustments() {
    let atom_each_chain = 10;
    let new_bond_length = 7.0;
    const n_chains: usize = 2;
    let n_atoms = atom_each_chain * n_chains;
    let mut system = SurpassAlphaSystem::new(&[atom_each_chain; n_chains], 1000.0);
    // ---------- Initialize coordinates
    let mut r = vec![3.8; n_atoms];
    let planar: Vec<f64> = (0..n_atoms).map(|_| 120.0_f64.to_radians()).collect();
    let dihedral: Vec<f64> = (0..n_atoms).map(|_| 180.0_f64.to_radians()).collect();
    let mut coords = vec![Vec3::default(); n_atoms];
    for i in 1..n_chains {
        r[i* atom_each_chain] = 10.0;
    }
    restore_linear_chain(&r[0..n_atoms], &planar[0..n_atoms], &dihedral[0..n_atoms], &mut coords[0..n_atoms]);
    for i in 0..n_atoms {
        system.vec3_to_ca(i,&coords[i]);
    }
    system.to_pdb_file("staring.pdb", false);
    system.adjust_bond_length(new_bond_length);
    system.to_pdb_file("adjusted.pdb", false);
    for ichain in 0..n_chains {
        for ires in 2..atom_each_chain {
            let a = system.atom_to_vec3(ichain* atom_each_chain + ires-2);
            let b = system.atom_to_vec3(ichain* atom_each_chain + ires-1);
            let c = system.atom_to_vec3(ichain* atom_each_chain + ires);
            assert_delta!(120.0, planar_angle3(&a, &b, &c).to_degrees(), 0.0001);
            assert_delta!(new_bond_length, b.distance_to(&a), 0.0001);
        }
    }
}