use clap::Result;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use crate::pdb_atom::PdbAtom;
use crate::pdb_compound::PdbCompound;
use crate::pdb_header::PdbHeader;
use crate::pdb_sequence_of_residue::PdbSequenceOfResidue;
use crate::pdb_source::PdbSource;
use crate::pdb_title::PdbTitle;

pub struct Pdb {
    pub header: Option<PdbHeader>,
    pub title: Option<PdbTitle>,
    pub compound: Option<PdbCompound>,
    pub source: Option<PdbSource>,
    pub sequence_of_residue: Option<PdbSequenceOfResidue>,
    atoms_list: Vec<PdbAtom>,
}

impl Pdb {
    pub fn new() -> Self {
        Self {
            header: None,
            title: None,
            compound: None,
            source: None,
            sequence_of_residue: None,
            atoms_list: vec![],
        }
    }

    pub fn from_file(file_name: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let file = File::open(file_name)?;
        let reader = BufReader::new(file);

        let mut pdb_file = Pdb::new();

        for line in reader.lines() {
            let line = line?;
            let splitted_line: Vec<&str> = line.split_whitespace().collect();

            match splitted_line[0] {
                "HEADER" => {
                    let header = PdbHeader::new(&line);
                    pdb_file.header = Some(header);
                },
                "TITLE" => {
                    let title = PdbTitle::new(line);
                    pdb_file.title = Some(title);
                },
                "COMPND" => {},
                "Source_" => {},
                "SequenceOfResidue_" => {},
                "ATOM" => {
                    let mut atom = PdbAtom::new(line);
                    atom.protein_name = pdb_file.header.as_ref().map_or_else(|| Path::new(file_name).file_stem().unwrap().to_str().unwrap().to_string(), |h| h.protein_name.clone());
                    pdb_file.atoms_list.push(atom);
                },
                _ => {},
            };
        }

        Ok(pdb_file)
    }

    pub fn write_csv(&self, file_path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let mut file = File::create(file_path)?;
        writeln!(file, "{}", PdbAtom::header())?;

        for atom in &self.atoms_list {
            writeln!(file, "{}", atom.to_csv_string())?;
        }

        Ok(())
    }

    pub fn write_to_file(&self, file_path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let mut file = File::create(file_path)?;
        let mut pdb_string = String::new();

        for atom in &self.atoms_list {
            pdb_string.push_str(&atom.to_string());
        }

        writeln!(file, "{}", pdb_string)?;

        Ok(())
    }

    pub fn get_atoms_list(&self) -> Vec<PdbAtom> {
        self.atoms_list.clone()
    }
}
