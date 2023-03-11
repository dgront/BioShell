use std::collections::HashMap;
use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};

#[derive(Default, Clone, Debug)]
/// Amino acid / nucleic sequence.
///
pub struct Sequence {
    /// identifies this sequence
    id: String,
    /// a sequence is represented as a vector of u8 bytes
    seq: Vec<u8>,
}

impl Sequence {
    /// Create a new instance of a Sequence from Strings.
    /// # Example
    /// ```rust
    /// use bioshell_core::sequence::Sequence;
    ///
    /// let seq_id: String = String::from("2gb1");
    /// let sequence: String = String::from("MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE");
    /// let seq = Sequence::new(&seq_id, &sequence);
    ///
    /// assert_eq!("MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE", seq.to_string())
    /// ```
    pub fn new(id: &String, seq: &String) -> Self {
        Sequence {
            id: id.to_owned(),
            seq: seq.chars().map(|c| c as u8).collect()
        }
    }

    /// Create a new instance of a Sequence by consuming the given data
    /// # Example
    /// ```rust
    /// use bioshell_core::sequence::Sequence;
    ///
    /// let seq_id = String::from("2gb1");
    /// let sequence: Vec<u8> = "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE".as_bytes().to_vec();
    /// let seq = Sequence::from_attrs(seq_id, sequence);
    ///
    /// let expected = "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE";
    /// assert_eq!(expected, seq.to_string());
    /// ```
    pub fn from_attrs(id: String, seq: Vec<u8>) -> Self {
        Sequence {id, seq}
    }

    /// Return the reference of the id of this Sequence.
    pub fn id(&self) -> &str { self.id.as_ref() }

    /// Return the reference of the sequence itself
    pub fn seq(&self) -> &Vec<u8> { &self.seq }

    /// Return the length of this sequence
    pub fn len(&self) -> usize { self.seq.len() }

    /// Returns an amino acid character at a given position in this `Sequence`
    pub fn char(&self, pos:usize) -> char { self.seq[pos] as char }

    /// Creates a string representing this sequence
    pub fn to_string(&self) -> String { String::from_utf8(self.seq.clone()).unwrap() }
}

impl fmt::Display for Sequence {
    /// Creates a `String` representation of a `Sequence` - FASTA format
    /// # Examples
    ///
    /// Create a `Sequence` and turn it into a string
    ///
    /// ```rust
    /// use bioshell_core::sequence::Sequence;
    /// use std::fmt::Write;
    /// // create a Sequence object
    /// let seq_id = String::from("2gb1");
    /// let sequence = "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE";
    /// let seq = Sequence::from_attrs(seq_id, sequence.to_vec());
    ///
    /// let mut actual = String::new();
    /// // populate `actual` with a string representation of the sequence
    /// write!(actual, "{}", seq).unwrap();
    ///
    /// let expected = "> 2gb1\nMTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE\n";
    ///
    /// assert_eq!(actual, expected)
    /// ```
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        write!(f, "> {}\n{}\n", self.id(), self.to_string())
    }
}

/// Read sequences in FASTA format from a buffer
///
/// To read a text file, use the [`from_file()`](from_file()) function
pub fn from_fasta_reader<R: Read>(reader: &mut R) -> Vec<Sequence> where R: BufRead {

    let mut out:Vec<Sequence> = Vec::new();
    let mut lines:Vec<String> = Vec::new();
    let mut buffer = String::new();

    loop {
        buffer.clear();
        let cnt = reader.read_line(&mut buffer);
        match cnt {
            Ok(cnt) => { if cnt == 0 { break; } },
            Err(e) => panic!("Cannot read from a FASTA buffer, error occurred:  {:?}", e),
        };
        lines.push(buffer.to_owned());
    }

    let mut header: String = String::new();
    let mut seq: String = String::new();
    for l in lines {
        let line: String = l.trim().to_owned();
        if line.len() == 0 { continue }
        if line.as_bytes()[0] as char == '>' {                  // --- It's a header!
            if seq.len() > 0 {                                  // --- create a new sequence entry and clear the string buffers
                out.push(Sequence::new(&header, &seq));
            }
            header = line[1..].to_owned();
            seq = String::new();
        } else {                                                // --- It's sequence
            seq = format!("{}{}",seq, line);
        }
    }
    if seq.len() > 0 {
        out.push(Sequence::new(&header, &seq));   // --- create the last sequence
    }

    return out;
}

/// Read sequences in Stockholm format
///
/// For the detailed description of the format, see: https://sonnhammer.sbc.su.se/Stockholm.html
pub fn from_stockholm_reader<R: Read>(reader: &mut R) -> Vec<Sequence> where R: BufRead {

    let mut out:Vec<Sequence> = Vec::new();
    let mut lines:Vec<String> = Vec::new();
    let mut buffer = String::new();
    let mut data: HashMap<String, String> = HashMap::new();
    let mut keys: Vec<String> = vec![];

    loop {
        buffer.clear();
        let cnt = reader.read_line(&mut buffer);
        match cnt {
            Ok(cnt) => { if cnt == 0 { break; } },
            Err(e) => panic!("Cannot read from a Stockholm buffer, error occurred:  {:?}", e),
        };
        if ! buffer.starts_with('#') {
            let tokens: Vec<&str> = buffer.split_ascii_whitespace().collect::<Vec<&str>>();
            if tokens.len() != 2 { continue;}
            else {
                let seq_id: String = tokens[0].to_owned();
                let seq: String = tokens[1].to_owned();
                match data.get_mut(&seq_id) {
                    None => {
                        data.insert(seq_id.clone(), seq);
                        keys.push(seq_id);
                    }
                    Some(subseq) => {
                        subseq.push_str(&seq);
                    }
                };
            }
            lines.push(buffer.to_owned());
        }
    }
    for key in keys {
        let seq: &String = data.get(&key).unwrap();
        out.push(Sequence::new(&key, seq));
    }

    return out;
}

/// Reads sequences from a file in a format known to BioShell
///
/// This function opens a given file, creates a buffered reader and passes it to a format-specific
/// reader function, such as [`from_stockholm_reader()`](from_stockholm_reader()) or
/// [`from_fasta_reader()`](from_fasta_reader())
pub fn from_file(filename: &str, reader_func: fn(&mut BufReader<File>) -> Vec<Sequence>) -> Vec<Sequence> {

    // Open the file in read-only mode (ignoring errors).
    let file = match File::open(filename) {
        Ok(file) => file,
        Err(_) => {
            eprintln!("\nCan't open an input  file: {}\n", filename);
            std::process::exit(1);
        }
    };
    let mut reader = BufReader::new(file);
    return reader_func(&mut reader);
}

/// Defines how A3M data will be processed
pub enum A3mConversionMode {
    /// insertions marked with lowercase characters will be removed
    RemoveSmallCaps,
    /// insertions marked with lowercase characters will be expanded with gaps in all other sequences in the given set
    ExpandWithGaps,
}

/// Process sequences read from A3M file by either expanding or removing insertions.
///
/// Lower case letters are used to denote parts of a sequence that is not aligned to another sequence.
/// A typical alignment in the *A3M* format may look like:
/// <pre>
/// MTYKLILNGKTLKgeTTTEAVDAAT
/// HQYKLALNGKTLKaqTTTEAVDAAT
/// </pre>
/// This alignment can be converted to *upper-case only* version (as in the FASTA format) by two ways:
///  - by aligning the lower-case letters with gaps; this would result in the following alignment:
/// <pre>
/// MTYKLILNGKTLKGE--TTTEAVDAAT
/// HQYKLALNGKTLK--AQTTTEAVDAAT
/// </pre>
///  - by removing them in both sequences (if possible), which would lead to:
/// <pre>
/// MTYKLILNGKTLKTTTEAVDAAT
/// HQYKLALNGKTLKTTTEAVDAAT
/// </pre>
/// This method removes the lower case letters accordingly to the desired scenario, requested by `mode` flag,
/// from all the sequences found in the given `sequences` vector.
///
/// # Examples
/// ```rust
/// use bioshell_core::sequence::{Sequence, a3m_to_fasta, A3mConversionMode};
/// // create two Sequence objects
/// let mut sequences = vec![Sequence::from_attrs(String::from("seq1"), "MTYKLILNGKTLKgeTTTEAVDAAT".to_vec()),
///         Sequence::new(String::from_attrs("seq2"), "HQYKLALNGKTLKaqTTTEAVDAAT".to_vec())];
///
/// a3m_to_fasta(&mut sequences, A3mConversionMode.RemoveSmallCaps);
/// assert_eq!(sequences[0].to_string(), "MTYKLILNGKTLKTTTEAVDAAT");
/// ```
pub fn a3m_to_fasta(sequences: &mut Vec<Sequence>, mode: &A3mConversionMode) {
    match mode {
        A3mConversionMode::RemoveSmallCaps => {
            for i in 0..sequences.len() {
                let seq = &sequences[i];
                let mut seq_data: Vec<u8> = vec![];
                for c in &seq.seq {
                    let cu = *c as char;
                    if cu.is_uppercase() || cu.eq(&'-') || cu.eq(&'_') {
                        seq_data.push(*c);
                    }
                }
                sequences[i] = Sequence::from_attrs(seq.id.clone(), seq_data);
            }
        }
        A3mConversionMode::ExpandWithGaps => { todo!() }
    }
}

/// Removes gaps from a sequence alignment
///
/// Given a gapped sequence and a (multiple) sequence alignment, this function removes
/// the columns from the alignment where the reference sequence has a gap.
/// # Arguments
/// * `reference` - a reference sequence
/// * `sequences` - aligned sequences, e.g. a multiple sequence alignment
///
/// # Example
/// ```rust
/// use bioshell_core::sequence::{Sequence, from_stockholm_reader, remove_gaps_by_sequence};
/// let msa = "A0A6G1KYC8/4-30                 -------TK----QT---TW-EK--PA------
/// A0A6G1KYC8/54-83                -------TK----AT---AW-EM--PK-QM---
/// A0A2Z6MTB8/60-86                -------TR----QS---SW-EK--P-------
/// A0A2Z6MTB8/102-129              -------TQ----QS---TW-TI--PEE-----
/// H2XMA8/16-49                    -------TQ----RT---TW-QD--PR------
/// UPI000A31515F/122-163           -------TR----ES---AW-TK--PD------
/// UPI000A31515F/383-441           -------TL----ES---TW-EK--PQE-----
/// UPI00057AF938/507-540           -------TR----TT---TW-KH--PC------
/// UPI00138FB958/985-1021          -------TQ----QT---SW-LH--PVSQ----";
/// let mut reader = BufReader::new(alignment.as_bytes());
/// let mut sequences = from_stockholm_reader(&mut reader);
/// let ref_seq = sequences[8].clone();
/// remove_gaps_by_sequence(&ref_seq, &mut sequences);
/// assert_eq!("TKQTTWEKPA--", sequences[0].to_string());
/// assert_eq!("TQQTSWLHPVSQ", sequences[8].to_string());
/// ```
pub fn remove_gaps_by_sequence(reference: &Sequence, sequences: &mut Vec<Sequence>) {

    let n_seq: usize = sequences.len();
    // --- check if all the sequences are of the same length
    for i in 0..n_seq {
        if sequences[i].len() != reference.len() {
            panic!("The following sequence has different length that the reference: {}", sequences[i].to_string());
        }
    }
    // --- create the list of indexes of elements to be copied
    let mut pos: Vec<usize> = vec![];
    for i in 0..reference.len() {
        let c = reference.seq[i] as char;
        if c != '-' && c != '_' { pos.push(i); }
    }
    // --- copy
    for j in 0..n_seq {
        let id = &sequences[j].id;
        let old_aa: &Vec<u8> = sequences[j].seq();
        let mut new_aa: Vec<u8> = vec![];
        for pi in &pos {
            new_aa.push(old_aa[*pi]);
        }
        sequences[j] = Sequence::from_attrs(id.clone(), new_aa);
    }
}