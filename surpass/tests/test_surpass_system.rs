
use std::io::BufReader;
use std::string::String;
use rand::rngs::SmallRng;
use rand::SeedableRng;
use bioshell_pdb::load_pdb_reader;
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
    assert_delta!(model.int_to_real(model.cax[0]), -13.296, 0.001);
    assert_delta!(model.int_to_real(model.cay[0]), 0.028, 0.001);
    assert_delta!(model.int_to_real(model.caz[0]), 3.924, 0.001);
    assert_delta!(model.int_to_real(model.cax[10]), -0.651, 0.001);
    assert_delta!(model.int_to_real(model.cay[10]), -2.752, 0.001);
    assert_delta!(model.int_to_real(model.caz[10]), 2.466, 0.001);

    // --- compute the distance between two atoms
    assert_delta!(model.distance(0, 1), 3.812, 0.01);
}

/// Make a system of 4 chains and check whether atoms are correctly assigned to chains
#[test]
fn test_4_chains() {

    let mut rnd = SmallRng::from_entropy();
    let model = SurpassAlphaSystem::make_random(&[10, 10, 10, 10], 100.0, &mut rnd);
    for i in 0..4 {
        assert_eq!(model.chain_atoms(i).start, i*10);
        assert_eq!(model.chain_atoms(i).end, (i+1)*10);
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
    assert_eq!(model.chain_atoms(0).start, 0);
    assert_eq!(model.chain_atoms(0).end, 10);
    for j in 0..10 {
        assert_eq!(model.chain(j), 0);
    }
}

#[test]
fn test_coords_operations() {
    let mut model = SurpassAlphaSystem::new(&[3], 100.0);

    model.cax[2] = model.real_to_int(6.0);
    assert_delta!(model.int_to_real(model.cax[2]), 6.0, 0.00001);
    // --- set X coordinate beyond the box
    model.cax[2] = model.real_to_int(51.0);
    assert_delta!(model.int_to_real(model.cax[2]), -49.0, 0.00001);
    // --- set X coordinate inside the box but close to the negative end
    model.cax[2] = model.real_to_int(-49.0);
    assert_delta!(model.int_to_real(model.cax[2]), -49.0, 0.00001);
    // --- set X coordinate outside the box on the negative side
    model.cax[2] = model.real_to_int(-51.0);
    assert_delta!(model.int_to_real(model.cax[2]), 49.0, 0.00001);
}

#[test]
fn test_diatance_evaluation() {
    let mut model = SurpassAlphaSystem::new(&[3], 100.0);
    for i in 0..3 {
        model.cax[i] = model.real_to_int(0.0);
        model.cay[i] = model.real_to_int(0.0);
        model.caz[i] = model.real_to_int(0.0);
    }
    model.cax[2] = model.real_to_int(6.0);
    model.cax[1] = model.real_to_int(3.0);
    model.cay[1] = model.real_to_int(4.0);
    assert_delta!(model.distance(0,1), 5.0, 0.00001);
    assert_delta!(model.distance(2,1), 5.0, 0.00001);
}