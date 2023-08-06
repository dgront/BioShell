use std::collections::HashSet;
use clap::Result;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::convert::TryFrom;
use std::fmt;

use crate::pdb_atom::PdbAtom;
use crate::pdb_compound::PdbCompound;
use crate::pdb_header::PdbHeader;
use crate::pdb_parsing_error::ParseError;
use crate::pdb_sequence_of_residue::PdbSequenceOfResidue;
use crate::pdb_source::PdbSource;
use crate::pdb_title::PdbTitle;
use crate::pdb_atom_filters::{SameResidue, PdbAtomPredicate, PdbAtomPredicate2};

pub struct ResidueIndex {
    pub chain_id: String,
    pub res_seq: i32,
    pub i_code: String
}

impl fmt::Display for ResidueIndex {

    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result { write!(f, "{}:{}{}", self.chain_id,self.res_seq,self.i_code) }
}

impl TryFrom<&PdbAtom> for ResidueIndex {
    type Error = ();

    fn try_from(a: &PdbAtom) -> Result<Self, Self::Error> {
        Ok(ResidueIndex{ chain_id: a.chain_id.clone(), res_seq: a.res_seq, i_code: a.i_code.clone() })
    }
}

// todo: Replace defined_sequence with a struct from bioshell-seq
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

    /// Provides immutable access to atoms of this  [`Structure`](Structure)
    pub fn atoms(&self) -> &Vec<PdbAtom> { &self.atoms }

    /// Provides mutable access to atoms of this  [`Structure`](Structure)
    pub fn atoms_mut(&mut self) -> &mut Vec<PdbAtom> { &mut self.atoms }

    pub fn chain_ids(&self) -> Vec<String> {
        let uniq: HashSet<&String> = self.atoms.iter().map(|a| &a.chain_id).collect();
        Vec::from_iter(uniq.iter().map(|s| *s).cloned())
    }

    pub fn residue_ids(&self) -> Vec<ResidueIndex> {

        self.atoms.iter().map(|a| ResidueIndex::try_from(a).unwrap()).collect()
    }

    pub fn residue_ids_by_chain(&self, chain_id: &str) -> Vec<ResidueIndex> {

        self.atoms.iter()
            .filter(|a| a.chain_id==chain_id)
            .map(|a| ResidueIndex::try_from(a).unwrap()).collect()
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




