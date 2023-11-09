#[cfg(test)]
mod tests {
    use bioshell_pdb::PdbHeader;

    #[test]
    fn test_pdb_header_new() {
        let header_line = "HEADER    KINASE                                  18-MAY-98   16PK              ".to_string();
        let header = PdbHeader::new(&header_line);
        assert_eq!(header.classification, "KINASE".to_string());
        assert_eq!(header.dep_date, "18-MAY-98".to_string());
        assert_eq!(header.id_code, "16PK".to_string());
    }

    #[test]
    fn test_pdb_header_to_string() {
        let header_line = "HEADER    KINASE                                  18-MAY-98   16PK              ".to_string();
        let header = PdbHeader::new(&header_line);
        let header_str = header.to_string();
        assert_eq!(header_str, "HEADER    KINASE                                  18-MAY-98   16PK".to_string());
    }
}