use std::fmt::Write;
use bioshell_core::sequence::{Sequence, from_fasta_reader, a3m_to_fasta, A3mConversionMode, from_stockholm_reader};
use std::io::BufReader;

#[test]
fn create_sequence() {
    // ---------- Create a simple amino acid sequence by calling new()
    {
        let read_id: String = String::from("2gb1");
        let sequence: String = String::from("MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE");
        let seq = Sequence::new(&read_id, &sequence);

        assert_eq!("MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE", seq.to_string())
    }
    {
        let read_id = String::from("2gb1");
        let sequence = "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE".as_bytes().to_vec();
        let seq = Sequence::from_attrs(read_id, sequence);

        let expected = "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE";
        assert_eq!(expected, seq.to_string());

        // ---------- Test the Display trait
        let expected_out = "> 2gb1\nMTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE\n";
        let mut actual = String::new();
        write!(actual, "{}", seq).unwrap();
        assert_eq!(actual, expected_out)
    }
}

#[test]
fn read_fasta() {
    let fasta = "> 2gb1
    MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE
    > 2azaA
AQCEATIESNDAMQYDLKEMVVDKSCKQFTVHLKHVGKMAKSAMGHNWVLTKEADKEGVATDGMNAGLAQDYVKAGDT
RVIAHTKVIGGGESDSVTFDVSKLTPGEAYAYFCSFPGHWAMMKGTLKLSN";
    let mut reader = BufReader::new(fasta.as_bytes());
    let records = from_fasta_reader(&mut reader);
    assert_eq!(2, records.len());
    assert_eq!(records[0].len(), 56);
}


#[test]
fn read_stockholm() {
    let input = "#=GF PDB 2gb1
#=GS DE 2gb1 protein
2gb1A                   MTYKLILNGKTLKGETTTEAVDAAT
UniRef100_UPI0000D834FD HQYKLALNGKTLKGETTTEAVDAAT

2gb1A                   AEKVFKQYANDNGVDGEWTYDDATK
UniRef100_UPI0000D834FD AEKVFKQYANDNGVDGEWTYDDATK

2gb1A                   TFTVTE
UniRef100_UPI0000D834FD TFTVTE
";
    let mut reader = BufReader::new(input.as_bytes());
    let records = from_stockholm_reader(&mut reader);
    assert_eq!(2, records.len());

    assert_eq!(records[0].len(), 56);
    assert_eq!(records[1].len(), 56);
    assert_eq!(records[0].id(), "2gb1A");
    assert_eq!(records[1].id(), "UniRef100_UPI0000D834FD");
    let ss: String = records[1].to_string();
    assert_eq!(ss, "HQYKLALNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE");
}

#[test]
fn convert_a3m() {
    let mut seqs = vec![Sequence::from_attrs(String::from("s1"), "ABC".as_bytes().to_vec()),
                    Sequence::from_attrs(String::from("s2"), "AaBC".as_bytes().to_vec()),
                    Sequence::from_attrs(String::from("s3"), "ABbC".as_bytes().to_vec())];
    a3m_to_fasta(&mut seqs, &A3mConversionMode::RemoveSmallCaps);
    let expected = "ABC";
    for s in seqs {
        assert_eq!(s.to_string(), expected);
    }
}