#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::BufReader;
    use bioshell_seq::msa::{longest_sequence, medoid_sequence, MSA, StockholmMSA};
    use bioshell_seq::SequenceError;

    #[allow(non_upper_case_globals)]
    const medoid_test: &str = ">Smed
MKTAYIA--KQRQISFVK--SHFTRQ
>L1
MKTAYIANNKQRVISFVKGGSHLTRQ
>L2
MKTAYIASSKQRQIAFVK--SHFARQ
>L3
MKTAYIA--KQREISFVKGGSHFTRQ
>L4
MKTAYIAQQKQRQISFAK--SHFTRQ";

    #[test]
    fn test_sequence_selectors() -> Result<(), SequenceError> {
        let file = File::open("tests/test_files/4Fe-4S-example.sto")?;
        let mut sto_reader = BufReader::new(file);
        let msa = StockholmMSA::from_stockholm_reader(&mut sto_reader)?;

        let seq = longest_sequence(&msa)?;
        assert_eq!(seq.len(), 80);

        let seq = medoid_sequence(&msa)?;
        assert_eq!(seq.len(), 80);

        let mut sto_reader = BufReader::new(medoid_test.as_bytes());
        let msa = MSA::from_fasta_reader(&mut sto_reader)?;
        let seq = medoid_sequence(&msa)?;
        assert_eq!(seq.seq(), b"MKTAYIAKQRQISFVKSHFTRQ");

        Ok(())
    }
}