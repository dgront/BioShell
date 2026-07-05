#[cfg(test)]
mod tests {
    use bioshell_seq::sequence::{parse_sequence_id, SeqId, SeqIdList};

    #[test]
    fn test_seq_id_detection() {
        let test_cases: Vec<(&str, SeqId)> = vec![
            (" PZ062902.1 ", SeqId::GenBank("PZ062902.1".to_string())),
            (" EU958043 ", SeqId::GenBank("EU958043".to_string())),
            ("[taxid=9606]", SeqId::TaxId("9606".to_string())),
            ("taxid=9606", SeqId::TaxId("9606".to_string())),
            ("TaxID=9606", SeqId::TaxId("9606".to_string())),
        ];

        for (input, expected) in test_cases {
            let parsed = parse_sequence_id(input);

            assert_eq!(parsed[0], expected);
            assert_eq!(parsed[0].to_string(), expected.to_string());
        }
    }

    #[test]
    fn test_seqidlist_display() {
        let ids = vec![
            SeqId::SwissProt("Q9NQX5".to_string()),
            SeqId::RefSeq("XP_001234567.1".to_string()),
            SeqId::PDB("1HHP:A".to_string()),
        ];

        let mut ids = SeqIdList::from(ids);
        ids.sort();
        let header = ids.to_string();
        assert_eq!(header, "pdb|1HHP:A|sp|Q9NQX5|ref|XP_001234567.1");

        let ids = vec![
            SeqId::SwissProt("Q9NQX5".to_string()),
            SeqId::RefSeq("XP_001234567.1".to_string()),
            SeqId::Organism("Homo Sapiens".to_string()),
            SeqId::TaxId("9606".to_string()),
        ];
        let mut ids = SeqIdList::from(ids);

        ids.sort();
        let header = ids.to_string();
        assert_eq!(header, "sp|Q9NQX5|ref|XP_001234567.1|taxid=9606 [organism=Homo Sapiens]");
    }

    #[test]
    fn test_seqid() {
        let mut ids = parse_sequence_id(">UniRef100_P81928 RPII140-upstream gene protein n=2 Tax=Drosophila melanogaster TaxID=7227 ");
        ids.sort();
        assert_eq!(ids.len(), 2);
        assert_eq!(ids[0], SeqId::UniRef("UniRef100_P81928".to_string()));

        let mut ids = parse_sequence_id(">sp.Q6GZX3.002L_FRG3G Uncharacterized protein OS=Frog virus 3 (isolate Goorha) OX=654924");
        ids.sort();
        assert_eq!("sp|Q6GZX3|002L_FRG3G|taxid=654924", &ids.to_string());
        assert_eq!(ids.len(), 3);
        assert_eq!(ids[0], SeqId::SwissProt("Q6GZX3".to_string()));
        assert_eq!(ids[1], SeqId::UniProtEntry("002L_FRG3G".to_string()));
        assert_eq!(ids[2], SeqId::TaxId("654924".to_string()));

        let mut ids = parse_sequence_id(">gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]");
        ids.sort();
        assert_eq!(ids.len(), 3);
        assert_eq!(ids[0], SeqId::GenBank("AAD44166.1".to_string()));
        assert_eq!(ids[1], SeqId::NCBIGI("5524211".to_string()));
        assert_eq!(ids[2], SeqId::Organism("Elephas maximus maximus".to_string()));
    }

    #[test]
    fn test_gi_id() {
        let ids = parse_sequence_id(">gi|5524211|");
        assert_eq!(ids[0], SeqId::NCBIGI("5524211".to_string()));
    }

    #[test]
    fn test_filename() {
        let test_cases: Vec<(&str, &str)> = vec![
            ("ref|XP_001234567.1|", "ref_XP_001234567.1"),
            ("blabla", "blabla"),
            ("sp|Q9NQX5|", "sp_Q9NQX5"),
            ("pdb|1HHP:A|", "pdb_1HHP"),
            ("pdb_00001abc", "pdb_00001abc"),
            ("gi|5524211|", "gi_5524211"),
            ("gb|AAD44166.1|", "gb_AAD44166.1"),
            ("PZ062902.1 [organism=Arabidopsis thaliana]", "gb_PZ062902.1_Arabidopsis_thaliana"),
            ("UniRef100_P81928", "UniRef100_P81928"),
            ("sp.Q6GZX3.002L_FRG3G", "sp_Q6GZX3_002L_FRG3G"),
            ("[Arabidopsis thaliana]", "Arabidopsis_thaliana"),
            ("[organism=Arabidopsis thaliana]", "Arabidopsis_thaliana"),
        ];

        for (input, expected) in test_cases {
            let ids = parse_sequence_id(input);
            let fname = ids.file_name();
// eprintln!("input: {input}, fname: {fname}, expected: {expected}");
            assert_eq!(fname, expected, "unexpected file_name() for input: {input}");
        }
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
            ("pdb_00002gb1|", true, "pdb_00002gb1"),
            ("pdb_00002gb1 ", true, "pdb_00002gb1"),
            ("pdb_00002gb1:ABC|", true, "pdb_00002gb1:ABC"),
            ("|2gb1ABC_", false, ""),
            ("ABC2gb1 ", false, ""),
            (" 3defX ", false, ""),
            ("XYZ3def_", false, ""),
        ];
        for (text, is_match, expected) in test_cases.iter() {
            let ids = parse_sequence_id(text);
            assert_eq!(ids.len(), 1);
            if *is_match {
                assert!(matches!(&ids[0], SeqId::PDB(_expected)), "failed on test case {}", &text);
            } else {
                assert!(matches!(ids[0], SeqId::Default(_)), "failed on test case {}", &text);
            }
        }
    }

    #[test]
    fn test_species_detection() {
        let test_cases: Vec<(&str, bool, &str)> = vec![
            ("[Homo sapiens]", true, "Homo sapiens"),
            ("[organism=Homo sapiens]",true, "Homo sapiens"),
            ("[9606]", false, ""),
            ("[123abc]", false, ""),
        ];
        for (text, is_match, expected) in test_cases.iter() {
            let ids = parse_sequence_id(text);
            if *is_match {
                assert!(matches!(&ids[0], SeqId::Organism(_expected)), "failed on test case {}", &text);
            } else {
                assert!(matches!(&ids[0], SeqId::Default(_expected)), "failed on test case {}", &text);
            }
        }
    }
}