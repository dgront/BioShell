use std::collections::{HashMap, HashSet};
use clap::Result;
use std::convert::TryFrom;

use itertools::{Itertools};
use bioshell_seq::chemical::{ResidueType, ResidueTypeManager, ResidueTypeProperties, KNOWN_RESIDUE_TYPES};
use bioshell_seq::sequence::Sequence;

use crate::pdb_atom::PdbAtom;
use crate::pdb_header::PdbHeader;
use crate::pdb_parsing_error::ParseError;
use crate::pdb_title::PdbTitle;
use crate::pdb_atom_filters::{SameResidue, PdbAtomPredicate, PdbAtomPredicate2, SameChain, ByResidueRange};
use crate::ResidueId;
use crate::secondary_structure::SecondaryStructure;


/// A biomacromolecular structure composed of [`PdbAtom`](PdbAtom) objects.
///
/// A [`Structure`](Structure) struct holds all atoms in a `Vec<PdbAtom>` container; its implementation
/// provides methods to look at them in various ways.
///
/// # Creating a [`Structure`](Structure)
/// Typically one gets a [`Structure`](Structure) object by loading it from a file in the PDB format:
/// ```no_run
/// use bioshell_pdb::{load_pdb_file, Structure};
/// let strctr = load_pdb_file("2gb1.pdb").unwrap();
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
/// by [`pdb_atom_filters`](crate::pdb_atom_filters) module.
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
/// # Removing atoms, residues or chains
/// This can be easily done using the [`retain()`](std::vec::Vec::retain()) method of a [`Vec`](std::vec::Vec) struct combined with a respective
/// [`PdbAtomPredicate`](PdbAtomPredicate). An example below **removes water molecules**:
///
/// ```
/// # use bioshell_pdb::{PdbAtom, Structure};
/// use bioshell_pdb::pdb_atom_filters::{IsWater, PdbAtomPredicate};
/// # let mut strctr = Structure::new();
/// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"));
/// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    518  O   HOH A  69      25.155  27.554  29.987  1.00 21.91           O"));
/// let hoh = IsWater;
/// strctr.atoms_mut().retain(|a| !hoh.check(&a));
/// # assert_eq!(strctr.count_atoms(), 1);
/// ```
/// Another example removes all hydrogen atoms:
/// ```
/// # use bioshell_pdb::{PdbAtom, Structure};
/// use bioshell_pdb::pdb_atom_filters::{IsHydrogen, PdbAtomPredicate};
/// # let mut strctr = Structure::new();
/// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    149  CA  GLY A   9      10.920  -2.963   0.070  1.00  0.18           C"));
/// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    153  HA2 GLY A   9      10.848  -2.565  -0.927  1.00  0.20           H"));
/// let is_h = IsHydrogen;
/// strctr.atoms_mut().retain(|a| !is_h.check(&a));
/// # assert_eq!(strctr.count_atoms(), 1);
/// ```
///
pub struct Structure {
    pub header: Option<PdbHeader>,
    pub title: Option<PdbTitle>,
    pub(crate) ter_atoms: HashMap<String, ResidueId>,
    pub(crate) atoms: Vec<PdbAtom>,
}

impl Structure {
    /// Create a new empty [`Structure`] that contains no atoms.
    pub fn new() -> Self {
        Self {
            header: None,
            title: None,
            ter_atoms: Default::default(),
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

    /// Counts chains of this [`Structure`](Structure)
    /// ```
    /// # use bioshell_pdb::{PdbAtom, Structure};
    /// let mut strctr = Structure::new();
    /// strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA A  68      25.790  28.757  29.513  1.00 16.12           C"));
    /// strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA B  68      25.790  28.757  29.513  1.00 16.12           C"));
    ///
    /// assert_eq!(strctr.count_chains(), 2);
    /// ```
    pub fn count_chains(&self) -> usize {
        let same_chain = SameChain{};
        return self.atoms().windows(2).filter(|a| !same_chain.check(&a[0], &a[1])).count() + 1
    }

    /// Provides immutable access to an atom
    /// ```
    /// # use bioshell_pdb::{PdbAtom, Structure, ResidueId};
    /// # let pdb_lines = vec!["ATOM    514  N   ALA A  68      26.532  28.200  28.365  1.00 17.85           N",
    /// #                     "ATOM    515  CA  ALA A  68      25.790  28.757  29.513  1.00 16.12           C",
    /// #                     "ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N",
    /// #                     "ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"];
    /// # let atoms: Vec<PdbAtom> = pdb_lines.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
    /// # let strctr = Structure::from_iterator(atoms.iter());
    /// let a = strctr.atom(&ResidueId::new("A", 69, ' ')," CA ").unwrap();
    /// assert_eq!(a.name, " CA ");
    /// # assert_eq!(a.res_seq, 69);
    /// ```
    pub fn atom(&self, res_id: &ResidueId, name: &str) -> Option<&PdbAtom> {
        self.atoms.iter().find(|&a| res_id.check(a) && a.name == name)
    }

    /// Provides mutable access to an atom
    pub fn atom_mut(&mut self, res_id: &ResidueId, name: &str) -> Option<&mut PdbAtom> {
        self.atoms.iter_mut().find(|a| res_id.check(a) && a.name == name)
    }

    /// Provides immutable access to atoms of this [`Structure`](Structure)
    pub fn atoms(&self) -> &Vec<PdbAtom> { &self.atoms }

    /// Provides mutable access to atoms of this  [`Structure`](Structure)
    pub fn atoms_mut(&mut self) -> &mut Vec<PdbAtom> { &mut self.atoms }

    /// Creates a vector that holds string identifiers for all chains of this [`Structure`](Structure)
    ///
    /// The vector is sorted alphabetically, regardless the order of chains in this [`Structure`](Structure)
    pub fn chain_ids(&self) -> Vec<String> {
        let uniq: HashSet<&String> = self.atoms.iter().map(|a| &a.chain_id).collect();
        let mut output = Vec::from_iter(uniq.iter().map(|s| *s).cloned());
        output.sort();

        return output;
    }

    /// Returns atoms of a given chain
    ///
    /// ```
    /// # use bioshell_pdb::{PdbAtom, Structure};
    /// # let pdb_lines = vec!["ATOM    514  N   ALA A  68      26.532  28.200  28.365  1.00 17.85           N",
    /// #                     "ATOM    515  CA  ALA A  68      25.790  28.757  29.513  1.00 16.12           C",
    /// #                     "ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N",
    /// #                     "ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"];
    /// # let atoms: Vec<PdbAtom> = pdb_lines.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
    /// # let strctr = Structure::from_iterator(atoms.iter());
    /// let chain_A_atoms = strctr.atoms_in_chain("A");
    /// # assert_eq!(chain_A_atoms.len(),4);
    /// ```
    pub fn atoms_in_chain(&self, chain_id: & str) -> Vec<&PdbAtom> {
        self.atoms.iter().filter(move |&atm| atm.chain_id==chain_id).collect()
    }

    /// Provide an iterator over atoms from a given range of residues
    ///
    /// The atoms may belong to different chains.
    ///
    /// # Example
    /// ```
    /// # use std::io::BufReader;
    /// use bioshell_pdb::{load_pdb_reader, ResidueId};
    /// # #[allow(non_upper_case_globals)]
    /// const pdb_txt: &str =
    /// "ATOM      2  CA  MET A   1     -13.296   0.028   3.924  1.00  0.43           C
    /// ATOM     21  CA  THR A   2      -9.669  -0.447   4.998  1.00  0.19           C
    /// ATOM     35  CA  TYR A   3      -7.173  -2.314   2.811  1.00  0.08           C
    /// ATOM     56  CA  LYS A   4      -3.922  -3.881   4.044  1.00  0.10           C
    /// ATOM     78  CA  LEU A   5      -0.651  -2.752   2.466  1.00  0.11           C
    /// ATOM     97  CA  ILE A   6       2.338  -5.105   2.255  1.00  0.13           C
    /// ATOM      2  CA  MET B   1     -13.296   0.028   3.924  1.00  0.43           C
    /// ATOM     21  CA  THR B   2      -9.669  -0.447   4.998  1.00  0.19           C
    /// ATOM     35  CA  TYR B   3      -7.173  -2.314   2.811  1.00  0.08           C
    /// ATOM     56  CA  LYS B   4      -3.922  -3.881   4.044  1.00  0.10           C
    /// ATOM     78  CA  LEU B   5      -0.651  -2.752   2.466  1.00  0.11           C";
    /// let strctr = load_pdb_reader(BufReader::new(pdb_txt.as_bytes())).unwrap();
    /// let first = ResidueId::new("A", 4, ' ');
    /// let last = ResidueId::new("B", 2, ' ');
    /// let mut iterator = strctr.atom_in_range(first, last);
    /// assert_eq!(iterator.count(), 5);
    /// ```
    pub fn atom_in_range(&self, first_res: ResidueId, last_res: ResidueId) -> impl Iterator<Item = &PdbAtom> {
        let check = ByResidueRange::new(first_res, last_res);
        self.atoms.iter().filter(move |&a| check.check(a))
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
        // --- return empty vector if there are no atoms in the given iterator
        return Vec::new();
    }

    /// Creates a vector of [`ResidueId`](ResidueId) object for each residue of this [`Structure`](Structure)
    ///
    /// This method simply calls [`Structure::residue_ids_from_atoms()`](Structure::residue_ids_from_atoms()) for `self` atoms
    pub fn residue_ids(&self) -> Vec<ResidueId> {

        Structure::residue_ids_from_atoms(self.atoms.iter())
    }

    /// Returns the type of a residue.
    ///
    /// The type is provided as a [`ResidueType`] object.
    /// ```
    /// # use bioshell_pdb::{PdbAtom, ResidueId, Structure};
    /// use bioshell_seq::chemical::StandardResidueType;
    /// let pdb_lines = vec!["ATOM    515  CA  ALA A  68      25.790  28.757  29.513  1.00 16.12           C"];
    /// let atoms: Vec<PdbAtom> = pdb_lines.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
    /// let strctr = Structure::from_iterator(atoms.iter());
    /// let res_type = strctr.residue_type(&ResidueId::new("A", 68, ' ')).unwrap();
    /// assert_eq!(res_type.code3, "ALA");
    /// assert_eq!(res_type.parent_type, StandardResidueType::ALA);
    /// ```
    pub fn residue_type(&self, res_id: &ResidueId) -> Result<ResidueType, ParseError> {

        // --- check if such a residue has at least one atom in this struct
        if let Some(atom) = self.atoms().iter().find(|&a| res_id.check(a)) {
            if let Some(res_type) = KNOWN_RESIDUE_TYPES.lock().unwrap().by_code3(&atom.res_name) {
                return Ok(res_type.clone());    // --- atom exists and its residue type is known
            } else {                            // --- atom exists but its residue type is NOT known
                return Err(ParseError::UnknownResidueType { res_type: atom.res_name.clone()});
            }
        } else {                        // --- atom doesn't exist
            return Err(ParseError::NoSuchResidue { res_id: res_id.clone() });
        }
    }

    /// Provides the [`ResidueId`] of the last residue in a given chain
    ///
    /// Any chain may contain residues and atoms that are listed after the `TER` residue; these are
    /// not covalently connected to their chain and are considered ligands.
    ///
    /// If a chain does not provide a `TER` record, ID of its very last residue is returned
    pub fn ter_residue(&self, chain_id: &str) -> ResidueId {
        if let Some(res_id) = self.ter_atoms.get(chain_id) {
            return res_id.clone();
        } else {
            let last_at = self.atoms.iter().rfind(|&a| a.chain_id==chain_id).unwrap();
            return ResidueId::try_from(last_at).unwrap();
        }
    }

    /// Provides a sequence of a given chain.
    ///
    /// Residues listed after the `TER` record are not included in the sequence returned by this method
    pub fn sequence(&self, chain_id: &str) -> Sequence {
        let ter_resid = self.ter_residue(chain_id);
        let atms = self.residue_first_atoms(chain_id);
        let mut aa: Vec<u8> = vec![];
        for a in atms {
            if let Some(restype) = ResidueTypeManager::get().by_code3(&a.res_name) {
                aa.push(u8::try_from(restype.parent_type.code1()).unwrap());
            } else {
                aa.push(b'X');
            }
            if ter_resid.check(a) {break}
        }
        if let Some(header) = &self.header {
            return Sequence::from_attrs(format!("{}:{}", header.id_code, chain_id), aa);
        }
        return Sequence::from_attrs(format!(":{}", chain_id), aa);
    }

    /// Provides a secondary structure of a given chain.
    ///
    /// Residues listed after the `TER` record are not included in the sequence returned by this method
    pub fn secondary(&self, chain_id: &str) -> SecondaryStructure {
        let ter_resid = self.ter_residue(chain_id);
        let atms = self.residue_first_atoms(chain_id);
        let mut sec: Vec<u8> = vec![];
        for a in atms {
            sec.push(a.secondary_struct_type);
            if ter_resid.check(a) {break}
        }
        return SecondaryStructure::from_attrs(sec);
    }

    /// Creates a vector of [`ResidueId`](ResidueId) object for each residue of a given chain
    ///
    pub fn chain_residue_ids(&self, chain_id: &str) -> Vec<ResidueId> {

        Structure::residue_ids_from_atoms(self.atoms.iter().filter(|&a| a.chain_id==chain_id))
    }

    /// Returns atoms of a given residue
    ///
    /// ```
    /// # use bioshell_pdb::{PdbAtom, ResidueId, Structure};
    /// # let pdb_lines = vec!["ATOM    514  N   ALA A  68      26.532  28.200  28.365  1.00 17.85           N",
    /// #                     "ATOM    515  CA  ALA A  68      25.790  28.757  29.513  1.00 16.12           C",
    /// #                     "ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N",
    /// #                     "ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"];
    /// # let atoms: Vec<PdbAtom> = pdb_lines.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
    /// # let strctr = Structure::from_iterator(atoms.iter());
    /// let res_atoms = strctr.residue_atoms(&ResidueId::new("A", 68, ' '));
    /// # assert_eq!(res_atoms.iter().count(),2);
    /// ```
    pub fn residue_atoms(&self, residue_id: &ResidueId) -> Vec<&PdbAtom> {
        self.atoms.iter().filter(|&a| residue_id.check(a)).collect()
    }

    /// Returns atoms of a given residue
    ///
    /// ```
    /// # use bioshell_pdb::{load_pdb_reader, PdbAtom, ResidueId, Structure};
    /// # use std::io::BufReader;
    /// let pdb_lines = "ATOM   4378  CA  HIS D 146      14.229  -1.501  26.683  1.00 31.89           C
    /// TER    4388      HIS D 146
    /// HETATM 4562 FE   HEM D 148      -1.727   4.699  23.942  1.00 15.46          FE";
    ///
    /// let mut strctr = load_pdb_reader(BufReader::new(pdb_lines.as_bytes())).unwrap();
    /// strctr.drop_ligands();
    /// assert_eq!(strctr.count_atoms(), 1);
    /// ```
    pub fn drop_ligands(&mut self) {
        for chain_id in self.chain_ids() {
            // --- check if TER is set; otherwise we won't  drop anything
            if let Some(res_id) = self.ter_atoms.get(&chain_id) {
                // --- check if TER residue has any atoms; otherwise we won't  drop anything
                if let Some(last_ter_atom) = self.atoms.iter().rfind(|&a| res_id.check(a)) {
                    let start_idx = self.atoms.iter().position(|a| a == last_ter_atom).unwrap() + 1;
                    let last_chain_atom = self.atoms.iter().rfind(|&a| a.chain_id == chain_id).unwrap();
                    let stop_idx = self.atoms.iter().position(|a| a == last_chain_atom).unwrap() + 1;
                    self.atoms.drain(start_idx..stop_idx);
                }
            }
        }
    }

    /// Sorts all the atoms of this structure.
    pub fn sort(&mut self) { self.atoms.sort(); }

    /// Returns true if all the atoms of this structure are sorted
    pub fn is_sorted(&self) -> bool {
        self.atoms.windows(2).all(|w| w[0] <= w[1])
    }

    #[allow(dead_code)]
    /// Borrows the very last atom of this [`Structure`]
    fn last_atom(&self) -> &PdbAtom { &self.atoms[self.atoms.len()-1] }

    /// Borrows the first atom from each residue of this [`Structure`]
    fn residue_first_atoms(&self, chain_id: &str) -> Vec<&PdbAtom> {
        let same_res = SameResidue {};
        let mut ats: Vec<&PdbAtom> = self.atoms().windows(2)
                .filter(|a| !same_res.check(&a[0], &a[1]))
                .filter(|a| &a[0].chain_id==chain_id && &a[1].chain_id==chain_id)
                .map(|a| &a[1]).collect();
        ats.insert(0,self.atoms.iter().find(|&a| a.chain_id==chain_id).unwrap());

        return ats;
    }
}

#[test]
fn test_first_residue_atoms() {
    let lines: [&str; 6] = [
        "ATOM    515  CA  ALA A  68      25.790  28.757  29.513  1.00 16.12           C",
        "ATOM    518  CB  ALA A  68      25.155  27.554  29.987  1.00 21.91           C",
        "ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C",
        "ATOM    518  CB  ALA A  69      25.155  27.554  29.987  1.00 21.91           C",
        "ATOM    515  CA  ALA B  68      25.790  28.757  29.513  1.00 16.12           C",
        "ATOM    518  CB  ALA B  68      25.155  27.554  29.987  1.00 21.91           C"];
    let atoms: Vec<PdbAtom> = lines.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
    let strctr = Structure::from_iterator(atoms.iter());

    assert_eq!(strctr.residue_first_atoms("A").len(), 2);
    assert_eq!(strctr.residue_first_atoms("B").len(), 1);
}

