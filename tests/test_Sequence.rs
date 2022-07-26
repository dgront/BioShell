use bioshell_core::Sequence;
use std::fmt::Write;
use bioshell_core::sequence::from_fasta_reader;
use std::io::BufReader;

#[test]
fn create_sequence() {
    // ---------- Create a simple amino acid sequence by calling new()
    {
        let read_id: String = String::from("2gb1");
        let sequence: String = String::from("MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE");
        let seq = Sequence::new(&read_id, &sequence);

        assert_eq!("> 2gb1\nMTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE\n", seq.to_string())
    }
    {
        let read_id = "2gb1";
        let sequence = "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE";
        let seq = Sequence::from_attrs(read_id, sequence);

        let expected = "> 2gb1\nMTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE\n";
        assert_eq!(expected, seq.to_string());

        // ---------- Test the Display trait
        let mut actual = String::new();
        write!(actual, "{}", seq).unwrap();
        assert_eq!(actual, expected)
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
    assert_eq!(2, records.len())
}