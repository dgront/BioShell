pub struct PdbSheet {
    pub strand_1: i32,
    pub sheet_id_2: String,
    pub num_strands_3: i32,
    pub init_res_name_4: String,
    pub init_chain_id_5: char,
    pub init_seq_num_6: i32,
    pub init_insert_code_7: char,
    pub end_res_name_8: String,
    pub end_chain_id_9: char,
    pub end_seq_num_10: i32,
    pub end_insert_code_11: char,
    pub sense_12: i32,
}

impl PdbSheet {
    fn get_key(&self) -> String {
        format!(
            "{}:{}{}-{}{}:{}{}-{}{}:{}",
            self.strand_1,
            self.init_chain_id_5,
            self.init_seq_num_6,
            self.init_insert_code_7,
            self.end_chain_id_9,
            self.end_seq_num_10,
            self.end_insert_code_11,
            self.num_strands_3,
            self.sense_12,
            self.sheet_id_2
        )
    }
}