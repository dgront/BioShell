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
        eprintln!("{:?}", ids);
        assert_eq!(ids.len(), 2);
        assert_eq!(ids[0], SeqId::UniProtID("002L_FRG3G".to_string()));
        assert_eq!(ids[1], SeqId::TaxId("654924".to_string()));
    }

    #[test]
    fn test_pbd_id() {
        let test_cases: Vec<(&str, bool)> = vec![
            ("pdb2gb1 ", true),
            ("pd2gb1 ", false),
            (" 2gb1 ", true),
            ("|2gb1 ", true),
            ("|2gb1|", true),
            ("|2gb1A ", true),
            ("|2gb1ABC ", true),
            ("|2gb1:ABC ", true),
            ("|2gb1:ABC|", true),
            ("|2gb1ABC_", false),
            ("ABC2gb1 ", false),
            (" 3defX ", true),
            ("XYZ3def_", false),
        ];
        for (i, (text, expected)) in test_cases.iter().enumerate() {
            let ids = parse_sequence_id(text);
            assert_eq!(ids.len(), 1);
            if *expected {
                assert!(matches!(ids[0], SeqId::PDB(_)), "failed on test case {}", &text);
            } else {
                assert!(matches!(ids[0], SeqId::Default(_)), "failed on test case {}", &text);
            }
        }
    }
}