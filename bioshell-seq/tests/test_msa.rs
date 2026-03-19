#[cfg(test)]
mod tests {

    use std::fs::File;
    use std::io::BufReader;
    use bioshell_seq::msa::{GfEntryType, MSA, StockholmMSA};
    use bioshell_seq::sequence::Sequence;
    use bioshell_seq::SequenceError;

    #[allow(non_upper_case_globals)]
    static msa_stockholm: &'static str = "# STOCKHOLM 1.0
#=GF PDB 2gb1
#=GS DE 2gb1 protein
2gb1A                   MTYKLILNGKTLKGETTTEAVDAAT
UniRef100_UPI0000D834FD HQYKLALNGKTLKGETTTEAVDAAT

2gb1A                   AEKVFKQYANDNGVDGEWTYDDATK
UniRef100_UPI0000D834FD AEKVFKQYANDNGVDGEWTYDDATK

2gb1A                   TFTVTE
UniRef100_UPI0000D834FD TFTVTE
";

    #[allow(non_upper_case_globals)]
    static msa_fasta: &'static str = "> seq A
MTYKL

> seq B
MTYKI

> seq C
MFYK-";


    #[test]
    fn read_msa_fasta() {
        let mut fas_reader = BufReader::new(msa_fasta.as_bytes());
        let msa = MSA::from_fasta_reader(&mut fas_reader).unwrap();
        assert_eq!(msa.n_seq(), 3);
    }

    #[test]
    fn iterate_over_column() {
        let msa = MSA::from_sequences(vec![Sequence::from_str("seq-1", "PERF"),
                                           Sequence::from_str("seq-2", "PDRV"),
                                           Sequence::from_str("seq-3", "PERV")]).unwrap();
        let iter = msa.nth_column_iter(1);
        let col0: Vec<u8> = iter.collect();
        assert_eq!(col0, [b'E', b'D', b'E']);
    }

    #[test]
    fn msa_as_fasta() {
        let expected = "> seq-1
PERF
> seq-2
PERV
";

        let msa = MSA::from_sequences(vec![Sequence::from_str("seq-1", "PERF"),
                                           Sequence::from_str("seq-2", "PERV")]).unwrap();

        let try_fasta = format!("{}", msa);
        assert_eq!(expected, try_fasta);
    }

    #[test]
    fn read_stockholm_msa() -> Result<(), SequenceError> {
        let file = File::open("tests/test_files/4Fe-4S-example.sto")?;
        let mut sto_reader = BufReader::new(file);
        let msa = StockholmMSA::from_stockholm_reader(&mut sto_reader).unwrap();

        let ids: Vec<&str> = msa.sequences().iter().map(|s| s.description()).collect();
        assert_eq!(
            ids,
            vec![
                "syb:TZ53_18285",
                "646611275",
                "EGCR1_03845",
                "SPSE_0307",
                "BBR47_24240",
                "SporoP37_10355",
            ]
        );

        // --- GF annotations
        assert_eq!(
            msa.gf(GfEntryType::ID).expect("missing GF ID"),
            ["Ferredoxin_4Fe4S"]
        );
        assert_eq!(
            msa.gf(GfEntryType::AC).expect("missing GF AC"),
            ["FERREDOXIN_4FE4S"]
        );
        assert_eq!(
            msa.gf(GfEntryType::DE).expect("missing GF DE"),
            ["Ferredoxin family protein containing a putative 4Fe-4S cluster-binding motif"]
        );
        assert_eq!(
            msa.gf(GfEntryType::TP).expect("missing GF TP"),
            ["Family"]
        );
        assert_eq!(
            msa.gf(GfEntryType::AU).expect("missing GF AU"),
            ["BioShell v.4 example"]
        );

        println!("{:40}", msa);
        Ok(())
    }
}