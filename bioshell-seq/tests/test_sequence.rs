use std::fmt::Write;
use bioshell_seq::sequence::{Sequence, a3m_to_fasta, A3mConversionMode, remove_gaps_by_sequence, FastaIterator, StockholmIterator, SequenceProfile, ProfileColumnOrder};
use std::io::BufReader;
use std::iter::zip;
use bioshell_seq::msa::MSA;


#[test]
fn test_sequence() {
    let header = String::from("gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]");
    let sequence = b"LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLV";
    let seq = Sequence::from_attrs(header, sequence.to_vec());
    assert_eq!("gi|5524211|gb|AAD44166.1|", seq.id());
}

#[test]
fn create_sequence() {
    // ---------- Create a simple amino acid sequence by calling new()
    {
        let read_id: String = String::from("2gb1");
        let sequence: String = String::from("MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE");
        let seq = Sequence::new(&read_id, &sequence);

        assert_eq!("MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE", seq.to_string(0))
    }
    {
        let read_id = String::from("2gb1");
        let sequence = "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE".as_bytes().to_vec();
        let seq = Sequence::from_attrs(read_id, sequence);

        let expected = "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE";
        assert_eq!(expected, seq.to_string(0));

        // ---------- Test the Display trait
        let expected_out = "> 2gb1\nMTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE\n";
        let mut actual = String::new();
        write!(actual, "{}", seq).unwrap();
        assert_eq!(actual, expected_out)
    }
}

#[allow(non_upper_case_globals)]
static fasta: &'static str = "> 2gb1
    MTYKLILNGKTLKGETTTEAVDAATAEKVF
    KQYANDNGVDGEWTYDDATKTFTVTE

    > 2azaA
AQCEATIESN DAMQYDLKEM VVDKSCKQFT VHLKHVGKMA

KSAMGHNWVL TKEADKEGVA TDGMNAGLAQ DYVKAGDTRV
IAHTKVIGGG ESDSVTFDVS KLTPGEAYAY FCSFPGHWAM
MKGTLKLSN
>";

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

#[test]
fn iterate_fasta() {

    let reader = BufReader::new(fasta.as_bytes());
    let iter_fasta = FastaIterator::new(reader);

    let records: Vec<Sequence> = iter_fasta.into_iter().collect();

    assert_eq!(2, records.len());
    assert_eq!(records[0].len(), 56);
    assert_eq!(records[1].to_string(0), "AQCEATIESNDAMQYDLKEMVVDKSCKQFTVHLKHVGKMAKSAMGHNWVLTKEADKEGVATDGMNAGLAQDYVKAGDTRVIAHTKVIGGGESDSVTFDVSKLTPGEAYAYFCSFPGHWAMMKGTLKLSN");
}

#[test]
fn read_stockholm() {

    let mut reader = BufReader::new(stockholm.as_bytes());
    let iter_sto = StockholmIterator::new(&mut reader);
    let records: Vec<Sequence> = iter_sto.into_iter().collect();
    assert_eq!(2, records.len());

    assert_eq!(records[0].len(), 56);
    assert_eq!(records[1].len(), 56);
    assert_eq!(records[0].description(), "2gb1A");
    assert_eq!(records[1].description(), "UniRef100_UPI0000D834FD");
    let ss: String = records[1].to_string(0);
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
        assert_eq!(s.to_string(0), expected);
    }
}

#[test]
fn remove_gaps() {
    let alignment = "UniRef90_A0A6G1KYC8/4-30                 -------TK----QT---TW-EK--PA------
UniRef90_A0A6G1KYC8/54-83                -------TK----AT---AW-EM--PK-QM---
UniRef90_A0A2Z6MTB8/60-86                -------TR----QS---SW-EK--P-------
UniRef90_A0A2Z6MTB8/102-129              -------TQ----QS---TW-TI--PEE-----
UniRef90_H2XMA8/16-49                    -------TQ----RT---TW-QD--PR------
UniRef90_UPI000A31515F/122-163           -------TR----ES---AW-TK--PD------
UniRef90_UPI000A31515F/383-441           -------TL----ES---TW-EK--PQE-----
UniRef90_UPI00057AF938/507-540           -------TR----TT---TW-KH--PC------
UniRef90_UPI00138FB958/985-1021          -------TQ----QT---SW-LH--PVSQ----
UniRef90_Q8CGF7-2/122-163                -------TR----ES---AW-TK--PD------
UniRef90_Q8CGF7-2/385-443                -------TL----ES---TW-EK--PQE-----
UniRef90_A0A1I8AJG4/92-118               -------TK----QS---SW-TK--PD------
UniRef90_A0A1I8AJG4/124-150              -------TK----TT---TW-TL--PE------
UniRef90_M4E0Z8/202-231                  -------TK----QS---TW-EK--PVE-----
UniRef90_M4E0Z8/244-272                  -------TK----QS---TW-TM--PEE-----
UniRef90_UPI0010A0DD51/122-163           -------TR----ES---AW-TK--PD------
UniRef90_UPI0010A0DD51/389-447           -------TL----ES---TW-EK--PQE-----
UniRef90_UPI000F744328/122-163           -------TR----ES---AW-TK--PD------
UniRef90_UPI000F744328/389-447           -------TL----ES---TW-EK--PQE-----";
    let mut reader = BufReader::new(alignment.as_bytes());
    let mut sequences: Vec<Sequence> = StockholmIterator::new(&mut reader).into_iter().collect();
    assert_eq!(19, sequences.len());

    assert_eq!("-------TQ----QT---SW-LH--PVSQ----", sequences[8].to_string(0));
    let ref_seq = sequences[8].clone();
    remove_gaps_by_sequence(&ref_seq, &mut sequences);

    assert_eq!("TKQTTWEKPA--", sequences[0].to_string(0));
    assert_eq!("TQQTSWLHPVSQ", sequences[8].to_string(0));
}

#[test]
fn create_sequence_profile() {
    let seq = zip(vec!["s1","s2","s3","s4"].iter(),vec!["cttaga","ctcaga","ctaagg","cttaga"].iter())
        .map(|(d,s)| Sequence::from_str(*d,*s)).collect();
    let msa = MSA::from_sequences(seq).unwrap();
    let profile = SequenceProfile::new(ProfileColumnOrder::dna_standard(),&msa);
    assert_eq!(profile.fraction(0, 1), 1.0);
}

#[test]
fn sequence_to_string() {
    let mut sequence = Sequence::from_str("test_seq", "P-RF_");
    bioshell_seq::sequence::remove_gaps(&mut sequence);
    let l = sequence.len();
    assert_eq!(l, 3);
    assert_eq!(sequence.to_string(0), String::from("PRF"));
}

#[test]
fn test_sequence_to_string_multiline() {
    let sequence = Sequence::from_str("test_seq", "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE");
    let expected_output = "MTYKLILNGK\nTLKGETTTEA\nVDAATAEKVF\nKQYANDNGVD\nGEWTYDDATK\nTFTVTE";

    assert_eq!(sequence.to_string(10), expected_output);
}