use std::string::String;

pub(crate) struct AtomLineParser {
    pub serial: Vec<i32>,
    pub name: Vec<String>,
    pub alt_loc: Vec<String>,
    pub res_name: Vec<String>,
    pub chain_id: Vec<String>,
    pub res_seq: Vec<i32>,
    pub i_code: Vec<String>,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub z: Vec<f64>,
    pub occupancy: Vec<f64>,
    pub temp_factor: Vec<f64>,
    pub element: Vec<String>,
    pub charge: Vec<f64>,
    pub is_hetero: Vec<bool>
}

impl AtomLineParser {

    pub fn new() -> AtomLineParser {
        return AtomLineParser {        serial: vec![],
            name: vec![],
            alt_loc: vec![],
            res_name: vec![],
            chain_id: vec![],
            res_seq: vec![],
            i_code: vec![],
            x: vec![],
            y: vec![],
            z: vec![],
            occupancy: vec![],
            temp_factor: vec![],
            element: vec![],
            charge: vec![],
            is_hetero: vec![]}
    }

    pub fn parse(&mut self, pdb_line: &str, is_hetero: bool) {

        self.is_hetero.push(is_hetero);

        self.serial.push(pdb_line[6..11].trim().parse::<i32>().unwrap());
        self.name.push(pdb_line[12..16].to_string());
        self.alt_loc.push(pdb_line[16..17].to_string());
        self.res_name.push(pdb_line[17..20].to_string());
        self.chain_id.push(pdb_line[21..22].to_string());
        self.res_seq.push(pdb_line[22..26].trim().parse::<i32>().unwrap());
        self.i_code.push(pdb_line[26..27].to_string());
        self.x.push(pdb_line[30..38].trim().parse::<f64>().unwrap());
        self.y.push(pdb_line[38..46].trim().parse::<f64>().unwrap());
        self.z.push(pdb_line[46..54].trim().parse::<f64>().unwrap());
        self.occupancy.push(pdb_line[54..60].trim().parse::<f64>().unwrap());
        self.temp_factor.push(pdb_line[60..66].trim().parse::<f64>().unwrap());
        // let segment_id_symbol = &pdb_line[76..78];
        // let charge_on_the_atom = &pdb_line[78..80];

    }

    // pub fn assemble_atom(pdb_atom: &PdbAtom) -> String {
    //     let atom_serial_no = pdb_atom.atom_serial_no_1.unwrap_or(0);
    //     let residue_no = pdb_atom.residue_no_9.unwrap_or(0);
    //
    //     let occupancy = pdb_atom.occupancy_12.unwrap_or(0.0);
    //     let temperature_factor = pdb_atom.temperature_factor_13.unwrap_or(0.0);
    //
    //     let atom = format!("{:<6}{:>5} {:<4}{:<1}{:<3} {:<1}{:>4}{:>1}{:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}{:>10}{:<2}{:<2}",
    //                        "ATOM", atom_serial_no, pdb_atom.atom_symbol_2, pdb_atom.alt_loc_indicator_6, pdb_atom.residue_name_7,
    //                        pdb_atom.chain_name_8, residue_no, pdb_atom.insertion_code_10, pdb_atom.coordinate_11.x, pdb_atom.coordinate_11.y, pdb_atom.coordinate_11.z, occupancy,
    //                        temperature_factor, pdb_atom.segment_identifier_14, pdb_atom.segment_identifier_symbol_15, pdb_atom.charge_of_the_atom_16);
    //
    //     atom
    // }
}