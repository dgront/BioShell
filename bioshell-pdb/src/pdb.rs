use clap::Result;
use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::pdb_atom::PdbAtom;
use crate::pdb_compound::PdbCompound;
use crate::pdb_header::PdbHeader;
use crate::pdb_parsing_error::PdbParseError;
use crate::pdb_sequence_of_residue::PdbSequenceOfResidue;
use crate::pdb_source::PdbSource;
use crate::pdb_title::PdbTitle;

pub struct Pdb {
    pub header: Option<PdbHeader>,
    pub title: Option<PdbTitle>,
    pub compound: Option<PdbCompound>,
    pub source: Option<PdbSource>,
    pub sequence_of_residue: Option<PdbSequenceOfResidue>,
    pub atoms_list: Vec<PdbAtom>,
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
    pub fn from_file(file_name: &str) -> Result<Self, PdbParseError> {
        let file = File::open(file_name)?;
        let reader = BufReader::new(file);
        let mut pdb_file = Pdb::new();
        //let mut line_count = 0;
        for line in reader.lines() {
            let line = line?;
            //line_count = line_count + 1;
            //println!("{:?}, {:?}", line_count, line);
            // Check that the line has a valid PDB record type
            let record_type = &line[0..6];
            let record = record_type.trim();
            match record {
                "HEADER" => {
                    let header = PdbHeader::new(&line);
                    pdb_file.header = Some(header);
                },
                "TITLE" => {
                    let title = PdbTitle::new(line.as_str());
                    pdb_file.title = Some(title);
                },
                //"COMPND" => {},
                //"Source_" => {},
                //"SequenceOfResidue_" => {},
                "ATOM" => {
                    let mut atom = PdbAtom::parse(line.as_str());
                    atom.is_hetero_atom = false;
                    atom.protein_name = pdb_file.header.as_ref().unwrap().get_protein_name().to_string();
                    pdb_file.atoms_list.push(atom);
                },
                "HETATM" => {
                    let mut atom = PdbAtom::parse(line.as_str());
                    atom.is_hetero_atom = true;
                    atom.protein_name = pdb_file.header.as_ref().unwrap().get_protein_name().to_string();
                    pdb_file.atoms_list.push(atom);
                },
                _ => {},
            };
        }
        Ok(pdb_file)
    }

/*
    pub fn write_csv(&self, file_path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let mut file = File::create(file_path)?;
        writeln!(file, "{}", PdbAtom::header())?;

        for atom in &self.atoms_list {
            let mut fields = vec![
                atom.atom_serial_no.unwrap().to_string(),
                atom.atom_symbol.clone(),
                atom.alt_loc_indicator.to_string(),
                atom.residue_name.clone(),
                atom.chain_name.to_string(),
                atom.residue_no.unwrap().to_string(),
                atom.insertion_code.to_string(),
                atom.get_coordinate().x.to_string(),
                atom.get_coordinate().y.to_string(),
                atom.get_coordinate().z.to_string(),
                atom.occupancy.map(|f| f.to_string()).unwrap_or_default(),
                atom.temperature_factor.map(|f| f.to_string()).unwrap_or_default(),
                atom.segment_identifier_symbol.to_string(),
                atom.charge_of_the_atom.to_string(),
                atom.protein_name.clone(),
            ];

            // Remove empty fields from the end of the vector
            while fields.last().map(|s| s.is_empty()).unwrap_or(false) {
                fields.pop();
            }

            // Write the fields to the CSV file
            writeln!(file, "{}", fields.join(","))?;
        }

        Ok(())
    }
*/
    /*
    pub fn write_to_file(&self, file_path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let mut file = File::create(file_path)?;
        let mut pdb_string = String::new();

        for atom in &self.atoms_list {
            pdb_string.push_str(&atom.to_string());
        }

        writeln!(file, "{}", pdb_string)?;

        Ok(())
    }*/


    pub fn get_atoms_list(&self) -> Vec<PdbAtom> {

        return self.atoms_list.clone();
    }
}
