use crate::pdb_sheet::PdbSheet;

pub struct PdbSheetParser {}

impl PdbSheetParser {
    pub fn parse_line(line: &str) -> Option<PdbSheet> {
        if !line.starts_with("SHEET") {
            return None;
        }
        let strand_1 = line[7..10].trim().parse().ok()?;
        let sheet_id_2 = line[11..14].trim().to_string();
        let num_strands_3 = line[14..16].trim().parse().ok()?;
        let init_res_name_4 = line[17..20].trim().to_string();
        let init_chain_id_5 = line.chars().nth(21)?;
        let init_seq_num_6 = line[22..26].trim().parse().ok()?;
        let init_insert_code_7 = line.chars().nth(26)?;
        let end_res_name_8 = line[28..31].trim().to_string();
        let end_chain_id_9 = line.chars().nth(32)?;
        let end_seq_num_10 = line[33..37].trim().parse().ok()?;
        let end_insert_code_11 = line.chars().nth(37)?;
        let sense_12 = line[38..40].trim().parse().ok()?;
        Some(PdbSheet {
            strand_1,
            sheet_id_2,
            num_strands_3,
            init_res_name_4,
            init_chain_id_5,
            init_seq_num_6,
            init_insert_code_7,
            end_res_name_8,
            end_chain_id_9,
            end_seq_num_10,
            end_insert_code_11,
            sense_12,
        })
    }
}