use std::io::BufReader;
use bioshell_seq::msa::MSA;
use bioshell_seq::sequence::Sequence;

#[allow(non_upper_case_globals)]
static stockholm: &'static str ="#=GF PDB 2gb1
#=GS DE 2gb1 protein
2gb1A                   MTYKLILNGKTLKGETTTEAVDAAT
UniRef100_UPI0000D834FD HQYKLALNGKTLKGETTTEAVDAAT

2gb1A                   AEKVFKQYANDNGVDGEWTYDDATK
UniRef100_UPI0000D834FD AEKVFKQYANDNGVDGEWTYDDATK

2gb1A                   TFTVTE
UniRef100_UPI0000D834FD TFTVTE
";

#[allow(non_upper_case_globals)]
static fasta: &'static str = "> seq A
MTYKL

> seq B
MTYKI

> seq C
MFYK-";

#[test]
fn read_msa_stockholm() {
    let mut sto_reader = BufReader::new(stockholm.as_bytes());
    let msa = MSA::from_stockholm_reader(&mut sto_reader).unwrap();
    assert_eq!(msa.n_seq(), 2);
}

#[test]
fn read_msa_fasta() {
    let mut fas_reader = BufReader::new(fasta.as_bytes());
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