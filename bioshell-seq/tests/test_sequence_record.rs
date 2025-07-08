#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::BufReader;
    use std::path::Path;
    use bioshell_io::open_file;
    use bioshell_seq::sequence::{NCBIIterator, SwissProtFolderIterator, SwissProtIterator};

    #[test]
    fn test_SwissProtIterator() -> Result<(), std::io::Error> {
        let reader = open_file("tests/test_files/A0A068Q5V6.spt")?;
        let mut iterator = SwissProtIterator::new(reader);
        let record = iterator.next().expect("expected one record");
        assert_eq!(record.accession, "A0A068Q5V6");
        assert_eq!(record.id, "C7150_PRUMU");
        assert_eq!(record.full_name, "Cytochrome P450 71AU50");
        assert_eq!(record.taxid, Some(102107));
        assert_eq!(record.sequence, "MVWIWATIGLLALVHILQAWWKNKKKRLPPGPRGFPIFGSLHLLGEFPNKDLHRLARKYGDIMYMRLGLMPTIVISSPEAAELFLKTHDLVFASRPPHEGSKHISFGQKNLIFSEYGAYWRDTRKMCTIELLSNHKINSFKSMRREEVSLCVESIRAAANNRGVAVDLSDKVSSLSVDMSCRMVLGKKYRDEEFDERGFKSVVREAIQLASAPNLGDYIRFIAPLDLQGFTKRMKSVNKAFDNLFEKIIEEHLQPNDGERTMDFVDVMVGFMGSEESEYRIERPHIKAIMLDMLVASMDTSATTIEWALSELMRHPKAMKKVQKELENVVGLDKMVEESDLEKLDYLNMVVKETFRLHPVAPLLIPHASIEDCTVNGYHIPKKSRVLINVWAIGRDPNAWTDAEKFIPERFEGSSVDVRGNHFQLIPFGSGRRRCPGIQLGLTVVQLVLAQLVHCFDWELPNNMLPEELDMTEEFGLTVPRAKHLLAIPSYRLRKSA");

        let reader = open_file("tests/test_files/R4K4X3.spt")?;
        let mut iterator = SwissProtIterator::new(reader);
        let record = iterator.next().expect("expected one record");
        assert_eq!(record.accession, "R4K4X3");
        assert_eq!(record.id, "R4K4X3_CLOPA");
        assert_eq!(record.full_name, "Ferredoxin");
        assert_eq!(record.taxid, Some(86416));
        assert_eq!(record.sequence, "MKGFVDKDTCIGCGLCTSICPEVFIMDDKGKAERSKNEILETLVASAQEAATECPVNAITVE");

        let seq = record.sequence();
        assert_eq!(seq.description(),"tr.R4K4X3.R4K4X3_CLOPA Ferredoxin OS=Clostridium pasteurianum BC1. OX=86416");

        Ok(())
    }

    #[test]
    fn test_SwissProtFolderIterator() -> Result<(), std::io::Error> {
        let all_records = SwissProtFolderIterator::from_folder("tests/test_files/", "spt")?;

        for record in all_records {
            println!("{} ({}) [{} AA]", record.id, record.accession, record.sequence.len());
        }
        Ok(())
    }

    #[test]
    fn test_NCBIIterator() {
        let path = Path::new("tests/test_files/WP_015613896.gp");
        let file = File::open(path).expect("failed to open WP_015613896.gp test file");
        let reader = BufReader::new(file);

        let mut iterator = NCBIIterator::new(reader);
        let record = iterator.next().expect("expected one record");
        assert_eq!(record.accession, "WP_015613896");
        assert_eq!(record.id, "WP_015613896");
        assert_eq!(record.taxid, Some(1501));
        assert_eq!(record.sequence, "MKGFVDKDTCIGCGLCTSICPEVFIMDDKGKAERSKNEILETLVASAQEAATECPVNAITVE");
    }
}
