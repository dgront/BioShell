
use std::io::BufReader;
use std::string::String;
use bioshell_pdb::load_pdb_reader;
use surpass::SurpassAlphaSystem;

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

#[test]
fn system_from_pdb() {

    let strctr = load_pdb_reader(BufReader::new(pdb_txt.as_bytes())).unwrap();

    let model = SurpassAlphaSystem::from_pdb_structure(&strctr, 100.0);
    // --- check if chains are properly assigned
    assert_eq!(model.count_chains(), 2);
    assert_eq!(model.chain(7), 1);
    assert_eq!(*model.chain_id(1).unwrap(), String::from("B"));

    // --- check coordinates
    assert!((model.cax(0) + 13.296).abs() < 0.001);
    assert!((model.cay(0) - 0.028).abs() < 0.001);
    assert!((model.caz(0) - 3.924).abs() < 0.001);
    assert!((model.cax(10) + 0.651).abs() < 0.001);
    assert!((model.cay(10) + 2.752).abs() < 0.001);
    assert!((model.caz(10) - 2.466).abs() < 0.001);

    // --- compute the distance between two atoms
    println!("{}", model.distance(0, 1));
    assert!((model.distance(0, 1) - 3.812) < 0.01);
}

#[test]
fn build_new_system() {
    let model = SurpassAlphaSystem::new(&[10, 10, 10, 10], 100.0);
    model.to_pdb_file("s.pdb", false);
}
