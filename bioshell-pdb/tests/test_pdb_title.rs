#[cfg(test)]
mod tests {
    use bioshell_pdb::PdbTitle;

    #[test]
    fn test_pdb_title() {
        let line = "TITLE     PHOSPHOGLYCERATE KINASE FROM TRYPANOSOMA BRUCEI BISUBSTRATE           ";
        let expected_text = "PHOSPHOGLYCERATE KINASE FROM TRYPANOSOMA BRUCEI BISUBSTRATE";
        let pdb_title = PdbTitle::from_pdb_line(line);

        assert_eq!(pdb_title.to_string(), expected_text);
    }
}