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
    /// a sequence is represented as a vector of characters
    seq: Vec<char>,
}


impl Sequence {
    /// Create a new instance of a Sequence.
    /// # Example
    /// ```rust
    /// use bioshell_core::sequence::Sequence;
    ///
    /// let read_id: String = String::from("2gb1");
    /// let sequence: String = String::from("MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE");
    /// let seq = Sequence::new(&read_id, &sequence);
    ///
    /// assert_eq!("MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE", seq.to_string())
    /// ```
    pub fn new(id: &String, seq: &String) -> Self {
        Sequence {
            id: id.to_owned(),
            seq: seq.chars().collect()
        }
    }
    /// Create a new instance of a Sequence from `&str` content.
    /// # Example
    /// ```rust
    /// use bioshell_core::sequence::Sequence;
    ///
    /// let read_id = "2gb1";
    /// let sequence = "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE";
    /// let seq = Sequence::from_attrs(read_id, sequence);
    ///
    /// let expected = "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE";
    /// assert_eq!(expected, seq.to_string());
    /// ```
    pub fn from_attrs(id: &str, seq: &str) -> Self {
        Sequence {
            id: String::from(id),
            seq: seq.chars().collect()
        }
    }
    /// Return the reference of the id of this Sequence.
    pub fn id(&self) -> &str { self.id.as_ref() }

    /// Return the reference of the sequence itself
    pub fn seq(&self) -> &Vec<char> { &self.seq }

    /// Return the length of this sequence
    pub fn len(&self) -> usize { self.seq.len() }

    /// Returns amino acid character at a given position in this `Sequence`
    pub fn char(&self, pos:usize) -> char { self.seq[pos] }

    /// Creates a string representing this sequence
    pub fn to_string(&self) -> String { self.seq.iter().collect::<String>() }
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
    /// let seq_id = "2gb1";
    /// let sequence = "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE";
    /// let seq = Sequence::from_attrs(seq_id, sequence);
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

/// Process sequences read from A3M file by either expanding or removing insertions
///
pub fn a3m_to_fasta(sequences: &mut Vec<Sequence>, mode: &A3mConversionMode) {
    match mode {
        A3mConversionMode::RemoveSmallCaps => {
            for i in 0..sequences.len() {
                let seq = &sequences[i];
                let seq_str: String = seq.seq.iter().
                    filter(|c| c.is_uppercase() || (*c).eq(&'-') || (*c).eq(&'_')).collect();
                sequences[i] = Sequence::new(&seq.id, &seq_str);
            }
        }
        A3mConversionMode::ExpandWithGaps => { todo!() }
    }
}