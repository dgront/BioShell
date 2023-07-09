#[cfg(test)]
mod tests {
    use bioshell_pdb::pdb_title::PdbTitle;

    #[test]
    fn test_pdb_title() {
        let line = "TITLE     PHOSPHOGLYCERATE KINASE FROM TRYPANOSOMA BRUCEI BISUBSTRATE           ";
        let expected_text = "PHOSPHOGLYCERATE KINASE FROM TRYPANOSOMA BRUCEI BISUBSTRATE";
        let pdb_title = PdbTitle::new(line);

        assert_eq!(pdb_title.to_string(), expected_text);
    }
}