use std::cmp::Ordering;
use std::cmp::Ordering::Equal;
use std::convert::TryFrom;
use std::fmt;

use crate::pdb_atom::PdbAtom;
use crate::pdb_atom_filters::{PdbAtomPredicate};
use crate::PDBError;

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
/// # let strctr = Structure::from_atoms("1xyz", atoms);
/// let res_id = ResidueId::new("A", 68, ' ');
/// let get_res = ByResidue::new(res_id);
/// # let mut cnt = 0;
/// for atom in strctr.atoms().iter().filter(|a| get_res.check(&a)) {
///     // ... process atoms of the residue 68 of chain A
///     cnt += 1;
/// }
/// # assert_eq!(cnt, 2);
/// ```
#[derive(Clone, Debug, Eq, Hash, Ord)]
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

    /// Creates a new [`ResidueId`](ResidueId) from a [`PdbAtom`](PdbAtom)
    ///
    /// # Example
    /// ```
    /// # use bioshell_pdb::{PdbAtom, ResidueId};
    /// # fn main() -> Result<(), ()> {
    /// let atom = PdbAtom::from_atom_line("ATOM    391  CA  LEU A  51      12.088   9.803  13.653  1.00  9.53           C  ");
    /// let res_id = ResidueId::try_from(&atom)?;
    /// assert_eq!(res_id, ResidueId::new("A", 51, ' '));
    /// # Ok(())
    /// # }
    /// ```
    fn try_from(a: &PdbAtom) -> Result<Self, Self::Error> {
        Ok(ResidueId { chain_id: a.chain_id.clone(), res_seq: a.res_seq, i_code: a.i_code })
    }
}

impl TryFrom<&str> for ResidueId {
    type Error = PDBError;

    /// Creates a new [`ResidueId`](ResidueId) from a string definition
    ///
    /// # Example
    /// ```
    /// # use bioshell_pdb::{PdbAtom, PDBError, ResidueId};
    /// # fn main() -> Result<(), PDBError> {
    /// assert_eq!(ResidueId::try_from("A:-1")?, ResidueId::new("A", -1, ' '));
    /// assert_eq!(ResidueId::try_from("AB:123")?, ResidueId::new("AB", 123, ' '));
    /// assert_eq!(ResidueId::try_from("AA:-1C")?, ResidueId::new("AA", -1, 'C'));
    /// # Ok(())
    /// # }
    /// ```
    fn try_from(res_id: &str) -> Result<Self, Self::Error> {
        let mut parts = res_id.splitn(2, ':');
        let prefix = parts.next().ok_or(PDBError::ResidueIdParsingError{ residue_id: res_id.to_string() })?;
        let rest = parts.next().ok_or(PDBError::ResidueIdParsingError{ residue_id: res_id.to_string() })?;

        // Separate the number and the optional character
        let mut num_part = String::new();
        let mut chars = rest.chars().peekable();

        // Handle optional minus sign
        if let Some(&'-') = chars.peek() {
            num_part.push('-');
            chars.next();
        }

        while let Some(&c) = chars.peek() {
            if c.is_ascii_digit() {
                num_part.push(c);
                chars.next();
            } else {
                break;
            }
        }

        let number: i32 = num_part.parse().map_err(|_e| PDBError::ResidueIdParsingError{ residue_id: res_id.to_string() })?;
        let suffix = chars.next().unwrap_or(' ');   // If there's nothing left after digits, use ' ' as default

        Ok(ResidueId::new(prefix, number, suffix))
    }
}

impl PartialEq for ResidueId {
    /// Check whether two [`ResidueId`](ResidueId) objects are equal.
    ///
    /// Returns true when `seq_res`, `chain_id` and `i_code` of the two [`ResidueId`](ResidueId)s are identical
    fn eq(&self, other: &Self) -> bool {
        self.res_seq == other.res_seq && self.chain_id == other.chain_id && self.i_code == other.i_code
    }
}

impl PartialOrd for ResidueId {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.chain_id < other.chain_id { return Some(Ordering::Less) }
        if self.chain_id > other.chain_id { return Some(Ordering::Greater) }
        if self.res_seq < other.res_seq { return Some(Ordering::Less) }
        if self.res_seq > other.res_seq { return Some(Ordering::Greater) }
        if self.i_code < other.i_code { return Some(Ordering::Less) }
        if self.i_code > other.i_code { return Some(Ordering::Greater) }
        return Some(Equal);
    }
}

impl PdbAtomPredicate for ResidueId {
    fn check(&self, a: &PdbAtom) -> bool {
        a.chain_id == self.chain_id && a.res_seq == self.res_seq && a.i_code == self.i_code
    }
}
