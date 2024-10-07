#[cfg(test)]
mod tests {
    use std::io::BufReader;
    use bioshell_cif::read_cif_buffer;
    use bioshell_pdb::{load_cif_reader, PdbHelix, PdbSheet, SecondaryStructure};

    #[test]
    fn helices_from_cif() {
        let cif_data = include_str!("./test_files/2gb1.cif");
        let reader = BufReader::new(cif_data.as_bytes());
        let cif_data = read_cif_buffer(reader).unwrap();
        let helices = PdbHelix::from_cif_data(&cif_data[0]).unwrap();
        assert_eq!(helices.len(), 1);
        let helix = &helices[0];
        assert_eq!(helix.length, 15);
        assert_eq!(helix.ser_num, "1");
        assert_eq!(helix.init_chain_id, "A");
        assert_eq!(helix.end_chain_id, "A");
        assert_eq!(helix.init_seq_num, 22);
        assert_eq!(helix.end_seq_num, 36);
        assert_eq!(helix.init_i_code, ' ');
        assert_eq!(helix.end_i_code, ' ');

        let cif_data = include_str!("./test_files/2fdo.cif");
        let reader = BufReader::new(cif_data.as_bytes());
        let cif_data = read_cif_buffer(reader).unwrap();
        let helices = PdbHelix::from_cif_data(&cif_data[0]).unwrap();
        assert_eq!(helices.len(), 8);
    }

    #[test]
    fn sheets_from_cif() {
        let cif_data = include_str!("./test_files/2gb1.cif");
        let reader = BufReader::new(cif_data.as_bytes());
        let cif_data = read_cif_buffer(reader).unwrap();
        let sheets = PdbSheet::from_cif_data(&cif_data[0]).unwrap();
        assert_eq!(sheets.len(), 4);

        let cif_data = include_str!("./test_files/2fdo.cif");
        let reader = BufReader::new(cif_data.as_bytes());
        let cif_data = read_cif_buffer(reader).unwrap();
        let sheets = PdbSheet::from_cif_data(&cif_data[0]).unwrap();
        assert_eq!(sheets.len(), 12);
    }

    #[test]
    fn secondary_from_cif() {
        let cif_data = include_str!("./test_files/2gb1.cif");
        let reader = BufReader::new(cif_data.as_bytes());
        let strctr = load_cif_reader(reader);
        assert!(strctr.is_ok());
        let strctr = strctr.unwrap();

        assert_eq!(strctr.secondary("A").to_string(),
                   "CEEEEEECCCCCCEEEEEECCHHHHHHHHHHHHHHHCCCCCEEEEECCCCEEEEEC");
    }

    #[test]
    fn create_secondary_structure() {
        let sec_str = SecondaryStructure::new(10);
        assert_eq!(sec_str.len(), 10);
        assert_eq!(sec_str.hec_code(5), b'C');
    }
}