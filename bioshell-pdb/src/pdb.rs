use std::collections::HashSet;
use clap::Result;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::convert::TryFrom;
use std::fmt;

use itertools::{Itertools};

use crate::pdb_atom::PdbAtom;
use crate::pdb_compound::PdbCompound;
use crate::pdb_header::PdbHeader;
use crate::pdb_parsing_error::ParseError;
use crate::pdb_sequence_of_residue::PdbSequenceOfResidue;
use crate::pdb_source::PdbSource;
use crate::pdb_title::PdbTitle;
use crate::pdb_atom_filters::{SameResidue, PdbAtomPredicate, PdbAtomPredicate2};

/// Unique identifier for a residue
///
/// Such an ID may be used to access atom of a residue from a [`Structure`](Structure)
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
/// let res_id = ResidueId::new("A", 68, " ");
/// let get_res = ByResidue::new(res_id);
/// # let mut cnt = 0;
/// for atom in strctr.atoms().iter().filter(|a| get_res.check(&a)) {
///     // ... process atoms of the residue 68 of chain A
///     cnt += 1;
/// }
/// # assert_eq!(cnt, 2);
/// ```
pub struct ResidueId {
    pub chain_id: String,
    pub res_seq: i32,
    pub i_code: String
}

impl ResidueId {
    /// Creates a new [`ResidueId`](ResidueId) from its properties
    pub fn new(chain_id: &str, res_seq: i32, i_code: &str) -> ResidueId {
        ResidueId{
            chain_id: chain_id.to_string(),
            res_seq,
            i_code: i_code.to_string()
        }
    }
}

impl fmt::Display for ResidueId {

    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result { write!(f, "{}:{}{}", self.chain_id,self.res_seq,self.i_code) }
}

impl TryFrom<&PdbAtom> for ResidueId {
    type Error = ();

    fn try_from(a: &PdbAtom) -> Result<Self, Self::Error> {
        Ok(ResidueId { chain_id: a.chain_id.clone(), res_seq: a.res_seq, i_code: a.i_code.clone() })
    }
}

/// A biomacromolecular structure composed of [`PdbAtom`](PdbAtom) objects.
///
/// A [`Structure`](Structure) struct holds all atoms in a `Vec<PdbAtom>` container; its implementation
/// provides methods to look at them in various ways.
///
/// # Creating a [`Structure`](Structure)
/// Typically one gets a [`Structure`](Structure) object by loading it from a file in the PDB format:
/// ```no_run
/// use bioshell_pdb::{load_pdb, Structure};
/// let strctr = load_pdb("2gb1.pdb").unwrap();
/// ```
/// A [`Structure`](Structure) can be also created from an [`Iterator`](Iterator) over [`PdbAtom`](PdbAtom)s:
/// ```
/// use bioshell_pdb::{PdbAtom, Structure};
/// let pdb_lines = vec!["ATOM    514  N   ALA A  68      26.532  28.200  28.365  1.00 17.85           N",
///                      "ATOM    515  CA  ALA A  68      25.790  28.757  29.513  1.00 16.12           C",
///                      "ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N",
///                      "ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"];
/// let atoms: Vec<PdbAtom> = pdb_lines.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
/// let strctr = Structure::from_iterator(atoms.iter());
/// # assert_eq!(strctr.count_atoms(), 4);
/// ```
///
/// # Accessing its atoms
/// A [`Structure`](Structure) implements two methods that provide mutable  and immutable borrow
/// of the vector of its atoms: [`atoms()`](Structure::atoms()) and [`atoms_mut()`](Structure::atoms_mut()) respectively.
/// These can be easily filtered by [`PdbAtomPredicate`](PdbAtomPredicate) predicates provided
/// by [`pdb_atom_filters`](pdb_atom_filters) module.
/// ```
/// # use bioshell_pdb::{PdbAtom, Structure};
/// use bioshell_pdb::pdb_atom_filters::{IsCA, PdbAtomPredicate};
/// # let pdb_lines = vec!["ATOM    514  N   ALA A  68      26.532  28.200  28.365  1.00 17.85           N",
/// #                     "ATOM    515  CA  ALA A  68      25.790  28.757  29.513  1.00 16.12           C",
/// #                     "ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N",
/// #                     "ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"];
/// # let atoms: Vec<PdbAtom> = pdb_lines.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
/// # let strctr = Structure::from_iterator(atoms.iter());
/// let is_ca = IsCA;
/// # let mut n_ca = 0;
/// for ca in  strctr.atoms().iter().filter(|a| is_ca.check(&a)) {
///     // --- process each alpha carbon here
/// # n_ca += 1;
/// }
/// # assert_eq!(n_ca, 2);
/// ```
///
///
/// # Accessing its residues and chains
/// A  [`Structure`](Structure) struct hold only atoms; chains and residues are not stored explicitely.
/// A list of residue IDs can be created on-demand from atoms:
/// ```
/// # use bioshell_pdb::{PdbAtom, Structure};
/// use bioshell_pdb::pdb_atom_filters::{ByResidue, PdbAtomPredicate};
/// # let pdb_lines = vec!["ATOM    514  N   ALA A  68      26.532  28.200  28.365  1.00 17.85           N",
/// #                     "ATOM    515  CA  ALA A  68      25.790  28.757  29.513  1.00 16.12           C",
/// #                     "ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N",
/// #                     "ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"];
/// # let atoms: Vec<PdbAtom> = pdb_lines.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
/// # let strctr = Structure::from_iterator(atoms.iter());
/// let residue_ids = Structure::residue_ids_from_atoms(strctr.atoms().iter());
/// for res_id in residue_ids {
///     let res_filter = ByResidue::new(res_id);
///     # let mut cnt = 0;
///     for atom in strctr.atoms().iter().filter(|a| res_filter.check(&a)) {
///         // ... process atoms of a given residue one by one
///         # cnt += 1
///     }
///     # assert_eq!(cnt, 2);
/// }
/// ```
///
pub struct Structure {
    pub header: Option<PdbHeader>,
    pub title: Option<PdbTitle>,
    pub compound: Option<PdbCompound>,
    pub source: Option<PdbSource>,
    pub defined_sequence: Option<PdbSequenceOfResidue>,
    atoms: Vec<PdbAtom>,
}

impl Structure {
    pub fn new() -> Self {
        Self {
            header: None,
            title: None,
            compound: None,
            source: None,
            defined_sequence: None,
            atoms: vec![],
        }
    }

    /// Creates a new [`Structure`](Structure) by filling it with atoms from an iterator.
    ///
    /// Atoms provided by an iterator will be cloned.
    ///
    /// # Example
    ///```
    /// # use bioshell_pdb::{PdbAtom, Structure};
    /// # use bioshell_pdb::pdb_atom_filters::{IsBackbone, PdbAtomPredicate};
    /// let mut strctr = Structure::new();
    /// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N"));
    /// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"));
    /// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    516  C   ALA A  69      26.891  29.054  30.649  1.00 15.28           C"));
    /// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    517  O   ALA A  69      26.657  29.867  31.341  1.00 20.90           O"));
    /// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    518  CB  ALA A  69      25.155  27.554  29.987  1.00 21.91           C"));
    /// let bb = IsBackbone{};
    /// let bb_strctr = Structure::from_iterator(strctr.atoms().iter().filter(|a|bb.check(a)));
    /// # assert_eq!(bb_strctr.count_atoms(), 4);
    /// ```
    pub fn from_iterator<'a, T: Iterator+Clone>(iter: T) -> Structure
        where T: Iterator<Item=&'a PdbAtom> {

        let mut strctr = Structure::new();
        for a in iter { strctr.atoms.push(a.clone()) }

        return strctr;
    }

    /// Pushes a given [`PdbAtom`](PdbAtom) at the end of this [`Structure`](Structure)
    /// The given atom object will be consumed in the process
    pub fn push_atom(&mut self, a: PdbAtom) {
        self.atoms.push(a);
    }

    /// Counts atoms of this [`Structure`](Structure)
    /// ```
    /// # use bioshell_pdb::{PdbAtom, Structure};
    /// let mut strctr = Structure::new();
    /// strctr.push_atom(PdbAtom::from_atom_line("ATOM    514  N   ALA A  68      26.532  28.200  28.365  1.00 17.85           N"));
    /// strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA A  68      25.790  28.757  29.513  1.00 16.12           C"));
    /// assert_eq!(strctr.count_atoms(), 2);
    /// ```
    pub fn count_atoms(&self) -> usize { self.atoms.len() }

    /// Counts residues of this [`Structure`](Structure)
    /// ```
    /// # use bioshell_pdb::{PdbAtom, Structure};
    /// let mut strctr = Structure::new();
    /// strctr.push_atom(PdbAtom::from_atom_line("ATOM    514  N   ALA A  68      26.532  28.200  28.365  1.00 17.85           N"));
    /// strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA A  68      25.790  28.757  29.513  1.00 16.12           C"));
    /// strctr.push_atom(PdbAtom::from_atom_line("ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N"));
    /// strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"));
    ///
    /// assert_eq!(strctr.count_residues(), 2);
    /// ```
    pub fn count_residues(&self) -> usize {
        let same_res = SameResidue{};
        return self.atoms().windows(2).filter(|a| !same_res.check(&a[0], &a[1])).count() + 1
    }

    /// Provides immutable access to atoms of this [`Structure`](Structure)
    pub fn atoms(&self) -> &Vec<PdbAtom> { &self.atoms }

    /// Provides mutable access to atoms of this  [`Structure`](Structure)
    pub fn atoms_mut(&mut self) -> &mut Vec<PdbAtom> { &mut self.atoms }

    pub fn chain_ids(&self) -> Vec<String> {
        let uniq: HashSet<&String> = self.atoms.iter().map(|a| &a.chain_id).collect();
        Vec::from_iter(uniq.iter().map(|s| *s).cloned())
    }

    /// Creates a vector of [`ResidueId`](ResidueId) object for each residue found in a given vector of atoms
    ///
    /// ```
    /// # use bioshell_pdb::{PdbAtom, Structure};
    /// let pdb_lines = vec!["ATOM    514  N   ALA A  68      26.532  28.200  28.365  1.00 17.85           N",
    ///                      "ATOM    515  CA  ALA A  68      25.790  28.757  29.513  1.00 16.12           C",
    ///                      "ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N",
    ///                      "ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"];
    /// let atoms: Vec<PdbAtom> = pdb_lines.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
    /// let res_idx = Structure::residue_ids_from_atoms(atoms.iter());
    /// assert_eq!(format!("{}", res_idx[0]), "A:68 ");
    /// assert_eq!(format!("{}", res_idx[1]), "A:69 ");
    /// assert_eq!(res_idx.len(), 2);
    /// ```
    pub fn residue_ids_from_atoms<'a>(atoms: impl Iterator<Item = &'a PdbAtom>) -> Vec<ResidueId> {
        // --- predicate used to check whether we are entering a new residue
        let same_res = SameResidue{};
        // --- turn a given iterator into a peekable one
        let mut peek_iter = atoms.peekable();
        // --- take the first PdbAtom, if there is any
        let maybe_first = peek_iter.peek();

        if let Some(&first) = maybe_first {
            // --- take the ResidueIdx of the first PdbAtom
            let first_idx = ResidueId::try_from(first).unwrap();
            let mut ret: Vec<ResidueId> = peek_iter.tuple_windows()
                .filter(|(a,b)| !same_res.check(&a, &b))
                .map(|(_, b)| ResidueId::try_from(b).unwrap()).collect();

            // --- insert the first ResidueIdx which we actually juped over during the iteration above
            ret.insert(0, first_idx);

            return ret;
        }
        return Vec::new();
    }

    /// Creates a vector of [`ResidueId`](ResidueId) object for each residue of this [`Structure`](Structure)
    ///
    /// This method simply calls [`Structure::residue_ids_from_atoms()`](Structure::residue_ids_from_atoms()) for `self` atoms
    pub fn residue_ids(&self) -> Vec<ResidueId> {

        Structure::residue_ids_from_atoms(self.atoms.iter())
    }

    pub fn residue_ids_by_chain(&self, chain_id: &str) -> Vec<ResidueId> {

        self.atoms.iter()
            .filter(|a| a.chain_id==chain_id)
            .map(|a| ResidueId::try_from(a).unwrap()).collect()
    }

}

pub fn load_pdb(file_name: &str) -> Result<Structure,ParseError> {

    let file = File::open(file_name)?;
    let reader = BufReader::new(file);
    let mut pdb_structure = Structure::new();

    let mut atoms: Vec<PdbAtom> = vec![];

    for line in reader.lines() {
        let line = line?;

        // Check that the line has a valid PDB record type
        let record_type = &line[0..6];
        let record = record_type.trim();
        match record {
            "HEADER" => {
                let header = PdbHeader::new(&line);
                pdb_structure.header = Some(header);
            },
            "TITLE" => {
                let title = PdbTitle::new(line.as_str());
                pdb_structure.title = Some(title);
            },
            //"COMPND" => {},
            //"Source_" => {},
            //"SequenceOfResidue_" => {},
            "ATOM" => {
                atoms.push(PdbAtom::from_atom_line(&line));
            },
            "HETATM" => {
                atoms.push(PdbAtom::from_atom_line(&line));
            },
            _ => {},
        };
    }

    println!("{:} atoms loaded",atoms.len());

    pdb_structure.atoms = atoms;

    Ok(pdb_structure)
}




