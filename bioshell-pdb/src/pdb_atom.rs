use std::string::String;

pub struct PdbAtom {
    pub serial: i32,
    pub name: String,
    pub alt_loc: String,
    pub res_name: String,
    pub chain_id: String,
    pub res_seq: i32,
    pub i_code: String,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub occupancy: f64,
    pub temp_factor: f64,
    pub element: Option<String>,
    pub charge: Option<String>,
    pub is_hetero_atom: bool,
    pub secondary_struct_symbol: char,
}

impl PdbAtom {
    pub fn new() -> PdbAtom {
        PdbAtom {
            serial: 1,
            name: String::new(),
            alt_loc: String::new(),
            res_name: String::from("ALA"),
            chain_id: String::new(),
            res_seq: 1,
            i_code: String::new(),
            x: 0.0,
            y: 0.0,
            z: 0.0,
            occupancy: 0.0,
            temp_factor: 0.0,
            element: None,
            is_hetero_atom: false,
            secondary_struct_symbol: 'C',
            charge: None
        }
    }

    // pub fn to_string(&self) -> String{
    //     return PdbLineParser::assemble_atom(&self);
    // }

    pub fn from_atom_line(pdb_line: &str, is_hetero: bool) -> PdbAtom {

        let serial = pdb_line[6..11].trim().parse::<i32>().unwrap();
        let name = pdb_line[12..16].to_string();
        let alt_loc = pdb_line[16..17].to_string();
        let res_name = pdb_line[17..20].to_string();
        let chain_id = pdb_line[21..22].to_string();
        let res_seq = pdb_line[22..26].trim().parse::<i32>().unwrap();
        let i_code = pdb_line[26..27].to_string();
        let x = pdb_line[30..38].trim().parse::<f64>().unwrap();
        let y = pdb_line[38..46].trim().parse::<f64>().unwrap();
        let z = pdb_line[46..54].trim().parse::<f64>().unwrap();
        let occupancy = pdb_line[54..60].trim().parse::<f64>().unwrap();
        let temp_factor = pdb_line[60..66].trim().parse::<f64>().unwrap();
        return PdbAtom{
            serial,
            name,
            alt_loc,
            res_name,
            chain_id,
            res_seq,
            i_code,
            x,
            y,
            z,
            occupancy,
            temp_factor,
            element: None,
            charge: None,
            is_hetero_atom: is_hetero,
            secondary_struct_symbol: 'C'
        };
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