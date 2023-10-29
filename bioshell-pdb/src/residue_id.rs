use std::convert::TryFrom;
use std::fmt;

use crate::pdb_atom::PdbAtom;
use crate::pdb_atom_filters::{PdbAtomPredicate};

/// Unique identifier for a residue
///
/// Such an ID may be used to access atom of a residue from a [`Structure`](crate::Structure)
/// # Example
/// ```
/// use bioshell_pdb::ResidueId;
/// # use bioshell_pdb::{PdbAtom, Structure};
/// # use bioshell_pdb::pdb_atom_filters::{ByResidue, PdbAtomPredicate};
/// # let pdb_lines = vec!["ATOM    514  N   ALA A  68      26.532  28.200  28.365  1.00 17.85           N",
/// #                     "ATOM    515  CA  ALA A  68      25.790  28.757  29.513  1.00 16.12           C",
/// #                     "ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N",
/// #                     "ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"];
/// # let atoms: Vec<PdbAtom> = pdb_lines.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
/// # let strctr = Structure::from_iterator(atoms.iter());
/// let res_id = ResidueId::new("A", 68, ' ');
/// let get_res = ByResidue::new(res_id);
/// # let mut cnt = 0;
/// for atom in strctr.atoms().iter().filter(|a| get_res.check(&a)) {
///     // ... process atoms of the residue 68 of chain A
///     cnt += 1;
/// }
/// # assert_eq!(cnt, 2);
/// ```
#[derive(Clone, Debug)]
pub struct ResidueId {
    pub chain_id: String,
    pub res_seq: i32,
    pub i_code: char
}

impl ResidueId {
    /// Creates a new [`ResidueId`](ResidueId) from its properties
    pub fn new(chain_id: &str, res_seq: i32, i_code: char) -> ResidueId {
        ResidueId{
            chain_id: chain_id.to_string(),
            res_seq,
            i_code: i_code
        }
    }
}

impl fmt::Display for ResidueId {

    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result { write!(f, "{}:{}{}", self.chain_id,self.res_seq,self.i_code) }
}

impl TryFrom<&PdbAtom> for ResidueId {
    type Error = ();

    fn try_from(a: &PdbAtom) -> Result<Self, Self::Error> {
        Ok(ResidueId { chain_id: a.chain_id.clone(), res_seq: a.res_seq, i_code: a.i_code })
    }
}

impl PartialEq<Self> for ResidueId {
    /// Check whether two [`ResidueId`](ResidueId) objects are equal.
    ///
    /// Returns true when `seq_res`, `chain_id` and `i_code` of the two [`ResidueId`](ResidueId)s are identical
    fn eq(&self, other: &Self) -> bool {
        self.res_seq==other.res_seq && self.chain_id==other.chain_id && self.i_code == other.i_code
    }
}

impl Eq for ResidueId {

}
impl PdbAtomPredicate for ResidueId {
    fn check(&self, a: &PdbAtom) -> bool {
        a.chain_id == self.chain_id && a.res_seq == self.res_seq && a.i_code == self.i_code
    }
}

/// Extracts [`ResidueId`] object from a `TER` pdb line
///
/// # Examples
/// ```
/// use bioshell_pdb::residue_id_from_ter_record;
/// let id1 = residue_id_from_ter_record("TER    1187      LEU B  75 ");
/// # assert_eq!(id1.i_code, ' ');
/// let id2 = residue_id_from_ter_record("TER    1187      LEU B  75A");
/// # assert_eq!(id2.i_code, 'A');
/// let id3 = residue_id_from_ter_record("TER    1187      LEU B  75");
/// assert_eq!(format!("{}", id3),"B:75 ");
/// ```
pub fn residue_id_from_ter_record(ter_line: &str) -> ResidueId {
    let res_seq: i32 = ter_line[22..26].trim().parse().ok().unwrap();
    let icode = if ter_line.len() > 26 { ter_line[26..27].chars().next().unwrap() } else {' '};
    return ResidueId::new(ter_line[21..22].trim(), res_seq, icode);
}