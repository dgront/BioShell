use crate::PdbAtom;

/// Calculate the squared distance between two atoms.
pub fn distance_squared(ai: &PdbAtom, aj:&PdbAtom) -> f64 {
    return ai.pos.distance_square_to(&aj.pos);
}

/// Calculate the distance between two atoms.
///
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