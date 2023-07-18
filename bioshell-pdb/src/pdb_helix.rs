pub struct PdbHelix {
    pub serial_number_1: i32,
    pub helix_id_2: String,
    pub init_res_name_3: String,
    pub init_chain_id_4: char,
    pub init_seq_num_5: i32,
    pub init_insert_code_6: char,
    pub end_res_name_7: String,
    pub end_chain_id_8: char,
    pub end_seq_num_9: i32,
    pub end_insert_code_10: char,
    pub helix_type_11: i32,
    pub helix_length_12: i32,
}

impl PdbHelix {
    fn get_key(&self) -> String {
        format!(
            "{}:{}{}-{}{}:{}{}-{}{}:{}",
            self.serial_number_1,
            self.init_chain_id_4,
            self.init_seq_num_5,
            self.init_insert_code_6,
            self.end_chain_id_8,
            self.end_seq_num_9,
            self.end_insert_code_10,
            self.helix_type_11,
            self.helix_length_12,
            self.helix_id_2
        )
    }
}