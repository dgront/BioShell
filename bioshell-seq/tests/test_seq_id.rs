#[cfg(test)]
mod tests {
    use bioshell_seq::sequence::{parse_sequence_id, SeqId, SeqIdList};

    #[test]
    fn test_seqidlist() {
        let ids = vec![
            SeqId::SwissProt("Q9NQX5".to_string()),
            SeqId::RefSeq("XP_001234567.1".to_string()),
            SeqId::PDB("1HHP:A".to_string()),
        ];

        let mut ids = SeqIdList::from(ids);
        ids.sort();
        let header = ids.to_string();
        assert_eq!(header, "PDB|1HHP:A|sp|Q9NQX5|RefSeq|XP_001234567.1");
    }

    #[test]
    fn test_seqid() {
        let ids = parse_sequence_id(">UniRef100_P81928 RPII140-upstream gene protein n=2 Tax=Drosophila melanogaster TaxID=7227 ");
        assert_eq!(ids.len(), 2);
        assert_eq!(ids[0], SeqId::UniRef("UniRef100_P81928".to_string()));

        let ids = parse_sequence_id(">sp.Q6GZX3.002L_FRG3G Uncharacterized protein OS=Frog virus 3 (isolate Goorha) OX=654924");
        assert_eq!(ids.len(), 2);
        assert_eq!(ids[0], SeqId::UniProtID("002L_FRG3G".to_string()));
        assert_eq!(ids[1], SeqId::TaxId("654924".to_string()));

        let ids = parse_sequence_id(">gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]");
        assert_eq!(ids.len(), 2);
        assert_eq!(ids[0], SeqId::GenBank("AAD44166.1".to_string()));
        assert_eq!(ids[1], SeqId::NCBIGI("5524211".to_string()));
    }

    #[test]
    fn test_gi_id() {
        let ids = parse_sequence_id(">gi|5524211|");
        assert_eq!(ids[0], SeqId::NCBIGI("5524211".to_string()));
    }



    #[test]
    fn test_pbd_id() {
        let test_cases: Vec<(&str, bool, &str)> = vec![
            ("pdb2gb1 ", true, "2gb1"),
            ("pd2gb1 ", false, ""),
            (" 2gb1 ", true, "2gb1"),
            ("|2gb1 ", true, "2gb1"),
            ("|2gb1|", true, "2gb1"),
            ("|2gb1A ", false, ""),
            ("|5524211| ", false, ""),
            ("|2gb1ABC ", false, ""),
            ("|2gb1:ABC ", true, "2gb1:ABC"),
            ("|2gb1:ABC ", true, "2gb1:ABC"),
            ("|2gb1:ABC|", true, "2gb1:ABC"),
            ("|2gb1ABC_", false, ""),
            ("ABC2gb1 ", false, ""),
            (" 3defX ", false, ""),
            ("XYZ3def_", false, ""),
        ];
        for (text, is_match, expected) in test_cases.iter() {
            let ids = parse_sequence_id(text);
            assert_eq!(ids.len(), 1);
            if *is_match {
                assert!(matches!(&ids[0], SeqId::PDB(expected)), "failed on test case {}", &text);
            } else {
                assert!(matches!(ids[0], SeqId::Default(_)), "failed on test case {}", &text);
            }
        }
    }
}