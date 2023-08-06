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

/// todo: Implement adding new atoms into a structure; remember to sort them; maybe insert after binary search?
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

    pub fn from_iterator<'a, T: Iterator+Copy>(iter: &T) where T: Iterator<Item=&'a PdbAtom>{
        let mut strctr = Structure::new();
        strctr.atoms = iter.cloned().collect();

    }

    pub fn atoms(&self) -> &Vec<PdbAtom> { &self.atoms }

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




