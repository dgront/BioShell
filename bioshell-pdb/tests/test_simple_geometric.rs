use bioshell_pdb::{PdbAtom, ResidueId, Structure};
use bioshell_pdb::calc::{distance};

#[allow(non_upper_case_globals)]
const pdb_2gb1:  &str = include_str!("./test_files/2gb1.pdb");

#[test]
fn test_distances() {
    let lines: Vec<_> = pdb_2gb1.split("\n").filter(|&l|l.starts_with("ATOM")).collect();
    let atoms: Vec<PdbAtom> = lines.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
    let strctr = Structure::from_iterator("1xyz", atoms.iter());
    let ai = strctr.atom(&ResidueId::new("A", 14, ' '), " CA ").unwrap();
    let aj = strctr.atom(&ResidueId::new("A", 15, ' '), " CA ").unwrap();
    assert_eq!(ai.res_name, "GLY");
    assert_eq!(aj.res_name, "GLU");
    assert!((distance(ai, aj) - 3.8).abs() < 0.1);
}