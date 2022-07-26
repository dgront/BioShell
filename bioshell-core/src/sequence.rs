use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};

#[derive(Default, Clone, Debug)]
pub struct Sequence {
    id: String,
    seq: String,
}


/// Amino acid / nucleic sequence.
///
impl Sequence {
    /// Create a new instance of a Sequence.
    /// # Example
    /// ```rust
    /// use bioshell_core::Sequence;
    ///
    /// let read_id: String = String::from("2gb1");
    /// let sequence: String = String::from("MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE");
    /// let seq = Sequence::new(&read_id, &sequence);
    ///
    /// assert_eq!("> 2gb1\nMTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE\n", seq.to_string())
    /// ```
    pub fn new(id: &String, seq: &String) -> Self {
        Sequence {
            id: id.to_owned(),
            seq: seq.to_owned()
        }
    }
    /// Create a new instance of a Sequence from `&str` content.
    /// # Example
    /// ```rust
    /// use bioshell_core::Sequence;
    ///
    /// let read_id = "2gb1";
    /// let sequence = "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE";
    /// let seq = Sequence::from_attrs(read_id, sequence);
    ///
    /// let expected = "> 2gb1\nMTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE\n";
    /// assert_eq!(expected, seq.to_string());
    /// ```
    pub fn from_attrs(id: &str, seq: &str) -> Self {
        Sequence {
            id: String::from(id),
            seq: String::from(seq)
        }
    }
    /// Return the id of this Sequence.
    pub fn id(&self) -> &str { self.id.as_ref() }

    /// Return the sequence itself
    pub fn seq(&self) -> &str { self.seq.as_ref() }
}

impl fmt::Display for Sequence {
    /// Creates a `String` representation of a `Sequence` - FASTA format
    /// # Examples
    ///
    /// Create a `Sequence` and turn it into a string
    ///
    /// ```rust
    /// use bioshell_core::Sequence;
    /// use std::fmt::Write;
    /// // create a Sequence object
    /// let read_id = "2gb1";
    /// let sequence = "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE";
    /// let seq = Sequence::from_attrs(read_id, sequence);
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
        write!(f, "> {}\n{}\n", self.id(), self.seq())
    }
}

pub fn from_fasta_reader<R: Read>(reader: &mut R) -> Vec<Sequence> where R: BufRead {

    let mut out:Vec<Sequence> = Vec::new();
    let mut lines:Vec<String> = Vec::new();
    let mut buffer = String::new();
    // let mut cnt: Option<usize>;

    loop {
        buffer.clear();
        let cnt = reader.read_line(&mut buffer);
        match cnt {
            Ok(cnt) => { if cnt == 0 { break; } },
            Err(e) => panic!("cannot from a FASTA buffer, error occurred:  {:?}", e),
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

pub fn from_fasta_file(filename: &str) -> Vec<Sequence> {
    let mut out: Vec<Sequence> = Vec::new();

    // Open the file in read-only mode (ignoring errors).
    let file = File::open(filename).unwrap();
    let mut reader = BufReader::new(file);
    return from_fasta_reader(&mut reader);
}
