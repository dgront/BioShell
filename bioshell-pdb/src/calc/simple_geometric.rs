use crate::{PdbAtom, ResidueId, Structure};
use crate::calc::dihedral_angle4;

/// Calculate the squared distance between two atoms.
pub fn distance_squared(ai: &PdbAtom, aj:&PdbAtom) -> f64 {
    return ai.pos.distance_square_to(&aj.pos);
}

/// Calculate the distance between two atoms.
///
/// # Example
///```rust
/// use bioshell_pdb::PdbAtom;
/// use bioshell_pdb::calc::distance;
/// let ai = PdbAtom::from_atom_line("ATOM      2  CA  MET A   1     -13.296   0.028   3.924  1.00  0.43           C");
/// let aj = PdbAtom::from_atom_line("ATOM     21  CA  THR A   2      -9.669  -0.447   4.998  1.00  0.19           C");
/// let d = distance(&ai, &aj);
/// # assert!((d-3.8).abs()<0.1);
///```
pub fn distance(ai: &PdbAtom, aj:&PdbAtom) -> f64 {
    return ai.pos.distance_to(&aj.pos);
}

///
/// # Example
/// ```
/// use bioshell_pdb::{assert_delta, PdbAtom, ResidueId, Structure};
/// use bioshell_pdb::calc::phi;
/// let pdb_lines = ["ATOM    339  C   ASP A  22     -10.767   3.482  -3.375  1.00  0.20           C",
/// "ATOM    349  N   ALA A  23     -10.594   2.767  -4.455  1.00  0.19           N",
/// "ATOM    350  CA  ALA A  23      -9.237   2.226  -4.777  1.00  0.18           C",
/// "ATOM    351  C   ALA A  23      -8.287   3.342  -5.231  1.00  0.19           C"];
/// let atoms: Vec<PdbAtom> = pdb_lines.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
/// let strctr = Structure::from_iterator(atoms.iter());
/// let phi_val = phi(&strctr, &ResidueId::new("A", 23, ' '));
/// assert!(phi_val.is_some());
/// assert_delta!(phi_val.unwrap().to_degrees(), -71.925, 0.01);
/// ```
pub fn phi(strctr: &Structure, which_res: &ResidueId) -> Option<f64> {
    if let Some(range) =  strctr.atoms_for_residue.get(which_res) {
        if range.start==0 { return None; }
        let pos = strctr.residue_ids().binary_search(which_res);
        if pos.is_err() {return None;}
        let prev = &strctr.residue_ids()[pos.unwrap() - 1];
        let c_prev = strctr.atom(prev, " C  ");
        let n = strctr.atom(which_res, " N  ");
        let ca = strctr.atom(which_res, " CA ");
        let c = strctr.atom(which_res, " C  ");
        if c_prev.is_none() || n.is_none() || ca.is_none() || c.is_none() {
            return None;
        }
        return Some(dihedral_angle4(&c_prev.unwrap().pos, &n.unwrap().pos,
                                    &ca.unwrap().pos, &c.unwrap().pos));
    } else { return None; }
}

///
/// # Example
/// ```
/// use bioshell_pdb::{assert_delta, PdbAtom, ResidueId, Structure};
/// use bioshell_pdb::calc::psi;
/// let pdb_lines = ["ATOM    349  N   ALA A  23     -10.594   2.767  -4.455  1.00  0.19           N",
/// "ATOM    350  CA  ALA A  23      -9.237   2.226  -4.777  1.00  0.18           C",
/// "ATOM    351  C   ALA A  23      -8.287   3.342  -5.231  1.00  0.19           C",
/// "ATOM    359  N   ALA A  24      -8.828   4.348  -5.860  1.00  0.22           N "];
/// let atoms: Vec<PdbAtom> = pdb_lines.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
/// let strctr = Structure::from_iterator(atoms.iter());
/// let psi_val = psi(&strctr, &ResidueId::new("A", 23, ' '));
/// assert!(psi_val.is_some());
/// assert_delta!(psi_val.unwrap().to_degrees(), -31.158, 0.01);
/// ```
pub fn psi(strctr: &Structure, which_res: &ResidueId) -> Option<f64> {
    if let Some(range) =  strctr.atoms_for_residue.get(which_res) {
        if range.end == strctr.atoms().len() { return None; }
        let pos = strctr.residue_ids().binary_search(which_res);
        if pos.is_err() {return None;}
        let next = &strctr.residue_ids()[pos.unwrap() + 1];
        let n = strctr.atom(which_res, " N  ");
        let ca = strctr.atom(which_res, " CA ");
        let c = strctr.atom(which_res, " C  ");
        let n_next = strctr.atom(next, " N  ");
        if n_next.is_none() || n.is_none() || ca.is_none() || c.is_none() {
            return None;
        }
        return Some(dihedral_angle4(&n.unwrap().pos, &ca.unwrap().pos,
                                    &c.unwrap().pos, &n_next.unwrap().pos,));
    } else { return None; }
}