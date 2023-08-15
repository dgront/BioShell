use bioshell_pdb::{PdbAtom, Structure};

#[allow(non_upper_case_globals)]
const lines:  [&str;9] = [
    "ATOM    514  CA  MET A  60      26.532  28.200  28.365  1.00 17.85           N",
    "ATOM    515  CA  CYS A  61      25.790  28.757  29.513  1.00 16.12           C",
    "ATOM    516  CA  GLY A  62      26.891  29.054  30.649  1.00 15.28           C",
    "ATOM    517  CA  ILE A  63      26.657  29.867  31.341  1.00 20.90           O",
    "ATOM    518  CA  ALA A  64      25.155  27.554  29.987  1.00 21.91           C",
    "ATOM    514  CA  MET B  61      26.532  28.200  28.365  1.00 17.85           N",
    "ATOM    515  CA  ALA B  62      25.790  28.757  29.513  1.00 16.12           C",
    "ATOM    516  CA  CYS B  63      26.891  29.054  30.649  1.00 15.28           C",
    "ATOM    518  CA  ALW B  64      25.155  27.554  29.987  1.00 21.91           C"];

#[test]
fn test_sequence_from_structure() {
    let atoms: Vec<PdbAtom> = lines.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
    let strctr = Structure::from_iterator(atoms.iter());
    let seq = strctr.sequence("A");
    assert_eq!(seq.to_string(), "MCGIA");

    let seq = strctr.sequence("B");
    assert_eq!(seq.to_string(), "MACX");
}

#[allow(non_upper_case_globals)]
const pdb_2gb1:  &str = include_str!("./test_files/2gb1.pdb");

#[test]
fn test_2gb1_loading() {
    let lines_2gb1: Vec<_> = pdb_2gb1.split("\n").filter(|&l|l.starts_with("ATOM")).collect();
    let atoms: Vec<PdbAtom> = lines_2gb1.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
    let strctr = Structure::from_iterator(atoms.iter());
    assert_eq!(strctr.count_atoms(), 855);
    assert_eq!(strctr.count_residues(), 56);
    assert_eq!(strctr.count_chains(), 1);
}