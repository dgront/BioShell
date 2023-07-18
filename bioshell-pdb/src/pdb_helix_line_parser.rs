use crate::pdb_helix::PdbHelix;

pub struct PdbHelixParser {}

impl PdbHelixParser {
    pub fn parse_line(line: &str) -> Option<PdbHelix> {
        if !line.starts_with("HELIX") {
            return None;
        }
        let serial_number_1 = line[11..15].trim().parse().ok()?;
        let helix_id_2 = line[15..19].trim().to_string();
        let init_res_name_3 = line[19..22].trim().to_string();
        let init_chain_id_4 = line.chars().nth(21)?;
        let init_seq_num_5 = line[22..26].trim().parse().ok()?;
        let init_insert_code_6 = line.chars().nth(26)?;
        let end_res_name_7 = line[31..34].trim().to_string();
        let end_chain_id_8 = line.chars().nth(33)?;
        let end_seq_num_9 = line[34..38].trim().parse().ok()?;
        let end_insert_code_10 = line.chars().nth(38)?;
        let helix_type_11 = line[39..41].trim().parse().ok()?;
        let helix_length_12 = line[71..76].trim().parse().ok()?;
        Some(PdbHelix {
            serial_number_1,
            helix_id_2,
            init_res_name_3,
            init_chain_id_4,
            init_seq_num_5,
            init_insert_code_6,
            end_res_name_7,
            end_chain_id_8,
            end_seq_num_9,
            end_insert_code_10,
            helix_type_11,
            helix_length_12,
        })
    }
}