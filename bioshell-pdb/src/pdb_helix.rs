pub struct PdbHelix {
    pub ser_num: i32,
    pub helix_id: String,
    pub init_res_name: String,
    pub init_chain_id: char,
    pub init_seq_num: i32,
    pub init_i_code: char,
    pub end_res_name: String,
    pub end_chain_id: char,
    pub end_seq_num: i32,
    pub end_i_code: char,
    pub helix_class: i32,
    pub length: i32,
}

impl PdbHelix {

    pub fn from_helix_line(line: &str) -> PdbHelix {

        let ser_num = line[11..15].trim().parse().unwrap();
        let helix_id = line[15..19].trim().to_string();
        let init_res_name = line[19..21].trim().to_string();
        let init_chain_id = line.chars().nth(21).unwrap();
        let init_seq_num = line[22..26].trim().parse().unwrap();
        let init_i_code = line.chars().nth(26).unwrap();
        let end_res_name = line[31..34].trim().to_string();
        let end_chain_id = line.chars().nth(33).unwrap();
        let end_seq_num = line[34..38].trim().parse().unwrap();
        let end_i_code = line.chars().nth(38).unwrap();
        let helix_class = line[39..41].trim().parse().unwrap();
        let length = line[71..76].trim().parse().unwrap();
        PdbHelix {
            ser_num: ser_num,
            helix_id: helix_id,
            init_res_name: init_res_name,
            init_chain_id: init_chain_id,
            init_seq_num: init_seq_num,
            init_i_code: init_i_code,
            end_res_name: end_res_name,
            end_chain_id: end_chain_id,
            end_seq_num: end_seq_num,
            end_i_code: end_i_code,
            helix_class: helix_class,
            length: length,
        }
    }
}