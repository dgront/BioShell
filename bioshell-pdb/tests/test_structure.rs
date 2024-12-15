use std::io::BufReader;
use itertools::assert_equal;
use bioshell_pdb::{Deposit, PdbAtom, ResidueId, Structure, PDBError};
#[allow(non_upper_case_globals)]
const lines_ca:  [&str;9] = [
    "ATOM    514  CA  MET A  60      26.532  28.200  28.365  1.00 17.85           N",
    "ATOM    515  CA  CYS A  61      25.790  28.757  29.513  1.00 16.12           C",
    "ATOM    516  CA  GLY A  62      26.891  29.054  30.649  1.00 15.28           C",
    "ATOM    517  CA  ILE A  63      26.657  29.867  31.341  1.00 20.90           O",
    "ATOM    518  CA  ALA A  64      25.155  27.554  29.987  1.00 21.91           C",
    "ATOM    514  CA  MET B  61      26.532  28.200  28.365  1.00 17.85           N",
    "ATOM    515  CA  ALA B  62      25.790  28.757  29.513  1.00 16.12           C",
    "ATOM    516  CA  CYS B  63      26.891  29.054  30.649  1.00 15.28           C",
    "ATOM    518  CA  ALW B  64      25.155  27.554  29.987  1.00 21.91           C"];

#[allow(non_upper_case_globals)]
const  lines_ala: [&str; 5] = [
    "ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N",
    "ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C",
    "ATOM    516  C   ALA A  69      26.891  29.054  30.649  1.00 15.28           C",
    "ATOM    517  O   ALA A  69      26.657  29.867  31.341  1.00 20.90           O",
    "ATOM    518  CB  ALA A  69      25.155  27.554  29.987  1.00 21.91           C",
];

#[test]
fn test_sequence_from_structure() {
    let atoms: Vec<PdbAtom> = lines_ca.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
    let strctr = Structure::from_iterator("1xyz", atoms.iter());
    let seq = strctr.sequence("A");
    assert_eq!(seq.to_string(80), "MCGIA");

    let seq = strctr.sequence("B");
    assert_eq!(seq.to_string(80), "MACX");
}

#[test]
fn struct_from_pushed_atoms() {
    let mut strctr = Structure::new("1xyz");
    for line in lines_ala {
        let atom = PdbAtom::from_atom_line(line);
        strctr.push_atom(atom);
    }
    let seq = strctr.sequence("A");
    assert_eq!(seq.to_string(80), "A");

    let ca_result = strctr.atom(&ResidueId::new("A", 69, ' '), " CA ");
    let ca = ca_result.unwrap();
    assert_eq!(ca.serial, 515);
}

#[allow(non_upper_case_globals)]
const pdb_2gb1:  &str = include_str!("./test_files/2gb1.pdb");

#[test]
fn test_2gb1_loading() {
    let lines_2gb1: Vec<_> = pdb_2gb1.split("\n").filter(|&l|l.starts_with("ATOM")).collect();
    let atoms: Vec<PdbAtom> = lines_2gb1.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
    let strctr = Structure::from_iterator("1xyz", atoms.iter());
    assert_eq!(strctr.count_atoms(), 855);
    assert_eq!(strctr.count_residues(), 56);
    assert_eq!(strctr.count_chains(), 1);
}

#[test]
fn test_secondary() -> Result<(), PDBError> {
    let deposit = Deposit::from_pdb_reader(BufReader::new(pdb_2gb1.as_bytes()))?;
    let strctr: Structure = deposit.structure();
    let mut ss_vec: Vec<u8> = vec![];
    for res in strctr.residue_ids() {
        ss_vec.push(strctr.residue_secondary(res).unwrap().hec_code());
    }
    let ss_str = String::from_utf8_lossy(&ss_vec).to_string();
    assert_eq!(&ss_str, &strctr.secondary("A").to_string());
    Ok(())
}

#[allow(non_upper_case_globals)]
const pdb_txt: &str =
    "ATOM      2  CA  MET A   1     -13.296   0.028   3.924  1.00  0.43           C
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
fn test_loading_from_reader() {
    let deposit = Deposit::from_pdb_reader(BufReader::new(pdb_txt.as_bytes())).unwrap();
    let strctr = deposit.structure();
    let seq = strctr.sequence("A");
    assert_eq!(seq.to_string(80), "MTYKLI");

    let seq = strctr.sequence("B");
    assert_eq!(seq.to_string(80), "MTYKL");
}

#[test]
fn test_atoms_by_range() {

    let deposit = Deposit::from_pdb_reader(BufReader::new(pdb_txt.as_bytes())).unwrap();
    let strctr = deposit.structure();
    let first = ResidueId::new("A", 4, ' ');
    let last = ResidueId::new("B", 2, ' ');
    let iterator = strctr.atoms_in_range(first, last);
    assert_eq!(iterator.count(), 5);
}

#[test]
fn test_atoms_by_residue() {

    let lines_2gb1: Vec<_> = pdb_2gb1.split("\n").filter(|&l|l.starts_with("ATOM")).collect();
    let atoms: Vec<PdbAtom> = lines_2gb1.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
    let strctr = Structure::from_iterator("1xyz", atoms.iter());

    assert_eq!(strctr.atoms_in_residue(&ResidueId::new("A", 1, ' ')).unwrap().count(), 19);
    assert_eq!(strctr.atoms_in_residue(&ResidueId::new("A", 56, ' ')).unwrap().count(), 16);
}
