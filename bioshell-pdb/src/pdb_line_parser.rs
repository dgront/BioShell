pub struct PdbLineParser;

impl PdbLineParser {
    pub fn parse_atom(pdb_line: &str) -> Option<[String; 16]> {
        let mut elements = [String::new(); 16];

        let record_name = &pdb_line[0..6];
        let atom_serial_no = &pdb_line[6..11];
        let atom_name = &pdb_line[12..16];
        let alt_loc_indicator = &pdb_line[16..17];
        let residue_name = &pdb_line[17..20];
        let chain_id = &pdb_line[21..22];
        let residue_seq_no = &pdb_line[22..26];
        let code_for_insertion_of_residues = &pdb_line[26..27];
        let x_coord = &pdb_line[30..38];
        let y_coord = &pdb_line[38..46];
        let z_coord = &pdb_line[46..54];
        let occupancy = &pdb_line[54..60];
        let temperature_factor = &pdb_line[60..66];
        let segment_id = &pdb_line[72..76];
        let segment_id_symbol = &pdb_line[76..78];
        let charge_on_the_atom = &pdb_line[78..80];

        elements[0] = record_name.to_string();
        elements[1] = atom_serial_no.to_string();
        elements[2] = atom_name.to_string();
        elements[3] = alt_loc_indicator.to_string();
        elements[4] = residue_name.to_string();
        elements[5] = chain_id.to_string();
        elements[6] = residue_seq_no.to_string();
        elements[7] = code_for_insertion_of_residues.to_string();
        elements[8] = x_coord.to_string();
        elements[9] = y_coord.to_string();
        elements[10] = z_coord.to_string();
        elements[11] = occupancy.to_string();
        elements[12] = temperature_factor.to_string();
        elements[13] = segment_id.to_string();
        elements[14] = segment_id_symbol.to_string();
        elements[15] = charge_on_the_atom.to_string();

        Some(elements)
    }

    pub fn assemble_atom(pdb_atom: &PdbAtom) -> String {
        let atom = format!("{:<6}{:>5} {:<4}{:<1}{:<3} {:<1}{:>4}{:>1}{:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}{:>10}{:<2}{:<2}",
                           "ATOM", pdb_atom.atom_serial_no, pdb_atom.atom_symbol, pdb_atom.alt_loc_indicator, pdb_atom.residue_name,
                           pdb_atom.chain_name, pdb_atom.residue_no, pdb_atom.insertion_code, pdb_atom.coordinate.x, pdb_atom.coordinate.y, pdb_atom.coordinate.z, pdb_atom.occupancy,
                           pdb_atom.temperature_factor, pdb_atom.segment_identifier, pdb_atom.segment_identifier_symbol, pdb_atom.charge_of_the_atom);

        atom
    }
}