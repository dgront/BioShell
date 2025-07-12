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
    }
}