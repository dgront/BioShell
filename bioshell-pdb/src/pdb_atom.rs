use std::string::String;
use bioshell_numerical::vec3::Vec3;
use crate::pdb_line_parser::PdbLineParser;

#[derive(Clone)]
pub struct PdbAtom {
    pub protein_name: String,
    pub atom_serial_no: Option<i32>,
    pub atom_symbol: String,
    pub atom_position: String,
    pub atom_no_in_the_branch: Option<i32>,
    pub connected_to_atom_no_in_the_branch: Option<i32>,
    pub alt_loc_indicator: String,
    pub residue_name: String,
    pub chain_name: String,
    pub residue_no: Option<i32>,
    pub insertion_code: String,
    pub coordinate: Vec3,
    pub occupancy: Option<f64>,
    pub temperature_factor: Option<f64>,
    pub segment_identifier: String,
    pub segment_identifier_symbol: String,
    pub charge_of_the_atom: String,
    pub is_hetero_atom: bool,
}

impl PdbAtom {
    pub fn new() -> PdbAtom {
        PdbAtom {
            protein_name: String::new(),
            atom_serial_no: Some(0),
            atom_symbol: String::new(),
            atom_position: String::new(),
            atom_no_in_the_branch: Some(0),
            connected_to_atom_no_in_the_branch: Some(0),
            alt_loc_indicator: String::new(),
            residue_name: String::new(),
            chain_name: String::new(),
            residue_no: Some(0),
            insertion_code: String::new(),
            coordinate: Vec3::zero(),
            occupancy: Some(0.0),
            temperature_factor: Some(0.0),
            segment_identifier: String::new(),
            segment_identifier_symbol: String::new(),
            charge_of_the_atom: String::new(),
        }
    }

    pub fn to_string(&self) -> String{
        return PdbLineParser::assemble_atom(&self);
    }

    pub fn parse(atom_line: &str) -> PdbAtom {
        let elements = PdbLineParser::parse_atom(atom_line).unwrap();

        let mut atom = PdbAtom::new();

        atom.atom_serial_no = Some(elements[1].trim().parse().unwrap());
        let atom_info = elements[2].trim();

        ///////////////////////////////////////
        //OE1
        //CA
        let atom_info_len = atom_info.len();
        if atom_info_len > 0 {
            atom.atom_symbol = atom_info[0..1].to_string();
        }
        if atom_info_len > 1 {
            atom.atom_position = atom_info[1..2].to_string();
        }
        if atom_info_len > 2 {
            let integer = atom_info[2..3].parse().unwrap_or(0);
            atom.atom_no_in_the_branch = Some(integer);
        }
        if atom_info_len > 3 {
            atom.connected_to_atom_no_in_the_branch = Some(atom_info[3..4].parse().unwrap());
        }

        ///////////////////////////////////////
        atom.alt_loc_indicator = elements[3].trim().to_string();
        atom.residue_name = elements[4].trim().to_string();
        atom.chain_name = elements[5].trim().to_string();
        atom.residue_no = Some(elements[6].trim().parse().unwrap());
        atom.insertion_code = elements[7].trim().to_string();
        atom.coordinate = Vec3::new(
            elements[8].trim().parse().unwrap(),
            elements[9].trim().parse().unwrap(),
            elements[10].trim().parse().unwrap(),
        );
        atom.occupancy = Some(elements[11].trim().parse().unwrap());
        atom.temperature_factor = Some(elements[12].trim().parse().unwrap());
        atom.segment_identifier =
            if elements[13].trim().is_empty() { String::new() } else { elements[13].trim().to_string() };
        atom.segment_identifier_symbol = elements[14].trim().to_string();
        atom.charge_of_the_atom =
            if elements[15].trim().is_empty() { String::new() } else { elements[15].trim().to_string() };

        atom
    }

    pub fn header() -> String {
        "AtomSerialNo,AtomSymbol,AtomPosition,AtomNoInTheBranch,ConnectedToAtomNoInTheBranch,AltLocIndicator,ResidueName,ChainName,ResidueNo,InsertionCode,Coordinate.X,Coordinate.Y,Coordinate.Z,Occupancy,TemperatureFactor,SegmentIdentifier,ChargeOfTheAtom".to_string()
    }

    pub fn to_csv_string(&self) -> String
    {
        let mut csv_string = String::new();
        csv_string.push_str(&self.atom_serial_no.map_or("#".to_string(), |x| x.to_string()));
        csv_string.push(',');
        csv_string.push_str(if self.atom_symbol.is_empty() { "#" } else { &self.atom_symbol });
        csv_string.push(',');
        csv_string.push_str(if self.atom_position.is_empty() { "#" } else { &self.atom_position });
        csv_string.push(',');
        csv_string.push_str(&self.atom_no_in_the_branch.map_or("#".to_string(), |x| x.to_string()));
        csv_string.push(',');
        csv_string.push_str(&self.connected_to_atom_no_in_the_branch.map_or("#".to_string(), |n| {
            if n == 0 {
                "#".to_string()
            } else {
                n.to_string()
            }
        }));
        csv_string.push(',');
        csv_string.push_str(if self.alt_loc_indicator.is_empty() { "#" } else { &self.alt_loc_indicator });//#
        csv_string.push(',');
        csv_string.push_str(&self.residue_name);
        csv_string.push(',');
        csv_string.push_str(&self.chain_name);
        csv_string.push(',');
        csv_string.push_str(&self.residue_no.map_or("#".to_string(), |x| x.to_string()));
        csv_string.push(',');
        csv_string.push_str(if self.insertion_code.is_empty() { "#" } else { &self.insertion_code });//#
        csv_string.push(',');
        csv_string.push_str(&self.coordinate.x.to_string());
        csv_string.push(',');
        csv_string.push_str(&self.coordinate.y.to_string());
        csv_string.push(',');
        csv_string.push_str(&self.coordinate.z.to_string());
        csv_string.push(',');
        csv_string.push_str(&self.occupancy.map_or("#".to_string(), |x| x.to_string()));
        csv_string.push(',');
        csv_string.push_str(&self.temperature_factor.map_or("#".to_string(), |x| x.to_string()));
        csv_string.push(',');
        csv_string.push_str(if self.segment_identifier.is_empty() { "#" } else { &self.segment_identifier });
        csv_string.push(',');
        csv_string.push_str(if self.segment_identifier_symbol.is_empty() { "#" } else { &self.segment_identifier_symbol });
        csv_string.push(',');
        csv_string.push_str(&self.charge_of_the_atom);
        return csv_string;
    }
        //
        // pub fn get_protein_name(&self) -> &str {
        //     &self.protein_name
        // }
        //
        // pub fn get_atom_serial_no(&self) -> Option<i32> {
        //     self.atom_serial_no
        // }
        //
        // pub fn get_atom_symbol(&self) -> &str {
        //     &self.atom_symbol
        // }
        //
        // pub fn get_atom_position(&self) -> &str {
        //     &self.atom_position
        // }
        //
        // pub fn get_atom_no_in_the_branch(&self) -> Option<i32> {
        //     self.atom_no_in_the_branch
        // }
        //
        // pub fn get_connected_to_atom_no_in_the_branch(&self) -> Option<i32> {
        //     self.connected_to_atom_no_in_the_branch
        // }
        //
        // pub fn get_alt_loc_indicator(&self) -> &str {
        //     &self.alt_loc_indicator
        // }
        //
        // pub fn get_residue_name(&self) -> &str {
        //     &self.residue_name
        // }
        //
        // pub fn get_chain_name(&self) -> &str {
        //     &self.chain_name
        // }
        //
        // pub fn get_residue_no(&self) -> Option<i32> {
        //     self.residue_no
        // }
        //
        // pub fn get_insertion_code(&self) -> &str {
        //     &self.insertion_code
        // }
        //
        pub fn get_coordinate(&self) -> &Vec3 {
             &self.coordinate
        }
        //
        // pub fn get_occupancy(&self) -> Option<f64> {
        //     self.occupancy
        // }
        //
        // pub fn get_temperature_factor(&self) -> Option<f64> {
        //     self.temperature_factor
        // }
        //
        // pub fn get_segment_identifier(&self) -> &str {
        //     &self.segment_identifier
        // }
        //
        // pub fn get_segment_identifier_symbol(&self) -> &str {
        //     &self.segment_identifier_symbol
        // }
        //
        // pub fn get_charge_of_the_atom(&self) -> &str {
        //     &self.charge_of_the_atom
        // }
    //}
}