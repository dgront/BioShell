
use std::collections::HashMap;
use std::fmt;
use std::io::{BufRead, BufReader};
use std::marker::PhantomData;
use std::hash::{Hash, Hasher};
use crate::errors::SequenceError;
use crate::msa::MSA;

#[derive(Default, Clone, Debug, PartialEq)]
/// Amino acid / nucleic sequence.
///
/// In `Rust` a char data type takes four bytes, which is not needed for biological alphabets. Therefore
/// the [`Sequence`](Sequence) struct stores an amino acid or a nucleic sequence as `Vec<u8>`.
pub struct Sequence {
    /// identifies this sequence
    description: String,
    /// a sequence is represented as a vector of u8 bytes
    seq: Vec<u8>,
}

impl Sequence {
    /// Create a new instance of a Sequence from Strings.
    /// # Example
    /// ```rust
    /// use bioshell_seq::sequence::Sequence;
    ///
    /// let seq_id: String = String::from("2gb1");
    /// let sequence: String = String::from("MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE");
    /// let seq = Sequence::new(&seq_id, &sequence);
    ///
    /// assert_eq!("MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE", seq.to_string())
    /// ```
    pub fn new(description: &String, seq: &String) -> Self {
        Sequence {
            description: description.to_owned(),
            seq: Self::filter_to_u8(seq)
        }
    }

    /// Create a new instance of a Sequence by consuming the given data
    ///
    /// # Example
    /// ```rust
    /// use bioshell_seq::sequence::Sequence;
    ///
    /// let seq_id = String::from("2gb1");
    /// let sequence: Vec<u8> = b"MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE".to_vec();
    /// let seq = Sequence::from_attrs(seq_id, sequence);
    ///
    /// let expected = "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE";
    /// assert_eq!(expected, seq.to_string());
    /// ```
    pub fn from_attrs(description: String, seq: Vec<u8>) -> Self {
        Sequence { description, seq}
    }

    /// A handy way to create a new Sequence from `str` data
    ///
    /// # Example
    /// ```rust
    /// use bioshell_seq::sequence::Sequence;
    ///
    /// let seq = Sequence::from_str("2gb1", "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE");
    ///
    /// let expected = "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE";
    /// assert_eq!(expected, seq.to_string());
    /// ```
    pub fn from_str(description: &str, seq: &str) -> Self {
        Self { description: String::from(description), seq: seq.as_bytes().to_vec()}
    }

    /// Return the description line of this Sequence.
    ///
    /// # Example
    /// ```rust
    /// use bioshell_seq::sequence::Sequence;
    ///
    /// let header = String::from("gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]");
    /// let sequence = b"LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLV";
    /// let seq = Sequence::from_attrs(header, sequence.to_vec());
    /// assert_eq!(seq.description(), "gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]");
    /// ```
    pub fn description(&self) -> &str { self.description.as_ref() }

    /// Return a string slice holding the ID of this sequence
    ///
    /// According to the NCBI standard, the FASTA header line should provide the unique identifier.
    /// The ID should be no longer than 25 characters and cannot contain any spaces; its fields should
    /// be separated by vertical bar (`|`) characters
    ///
    /// # Example
    /// ```rust
    /// use bioshell_seq::sequence::Sequence;
    ///
    /// let header = String::from("gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]");
    /// let sequence = b"LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLV";
    /// let seq = Sequence::from_attrs(header, sequence.to_vec());
    /// assert_eq!("gi|5524211|gb|AAD44166.1|", seq.id());
    /// ```
    pub fn id(&self) -> &str { self.description.split_whitespace().next().unwrap() }

    /// Return the reference of the sequence itself
    pub fn seq(&self) -> &Vec<u8> { &self.seq }

    /// Return the length of this sequence
    pub fn len(&self) -> usize { self.seq.len() }

    /// Returns a residue (e.g. an amino acid or a nucleotide) character at a given position in this `Sequence`
    pub fn char(&self, pos:usize) -> char { self.seq[pos] as char }

    /// Returns a residue (e.g. an amino acid or a nucleotide) u8 code at a given position in this `Sequence`
    pub fn u8(&self, pos:usize) -> u8 { self.seq[pos] }

    /// Creates a string representing this sequence.
    ///
    /// The returned string contains only the sequence itself without description.
    /// # Example
    /// ```rust
    /// use bioshell_seq::sequence::Sequence;
    ///
    /// let header = String::from("gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]");
    /// let sequence = b"LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLV";
    /// let seq = Sequence::from_attrs(header, sequence.to_vec());
    /// assert_eq!(seq.to_string(), "LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLV");
    /// ```
    pub fn to_string(&self) -> String { String::from_utf8(self.seq.clone()).unwrap() }

    fn filter_to_u8(seq: &String) -> Vec<u8> {
        seq.chars().filter(|c| c != &' ' && c != &'*').map(|c| c as u8).collect()
    }
}

impl fmt::Display for Sequence {
    /// Creates a `String` representation of a `Sequence` - FASTA format
    /// # Examples
    ///
    /// Create a `Sequence` and turn it into a string
    ///
    /// ```rust
    /// use bioshell_seq::sequence::Sequence;
    /// use std::fmt::Write;
    /// // create a Sequence object
    /// let seq_id = String::from("2gb1");
    /// let sequence = b"MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE";
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
        write!(f, "> {}\n{}\n", self.description(), self.to_string())
    }
}


impl Hash for Sequence {
    /// Computes a hash for a given sequence
    ///
    /// The hash depends only on the sequence itself, other data such as a header is irrelevant here.
    /// This allows detection of identical sequences even if their headers differs, e.g in the case
    /// of two identical chains of a PDB deposit
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.seq.hash(state);
    }
}

/// Iterator that provides sequences from a FASTA-formatted buffer.
///
/// This object iterates over a buffer without loading its whole content which allows processing
/// very large FASTA files.
pub struct FastaIterator<R> {
    reader: BufReader<R>,
    buffer: String,
    header: String,
    seq: String,
}

impl<R: BufRead> FastaIterator<R> {
    pub fn new(stream: R) -> Self {
        FastaIterator {
            reader: BufReader::new(stream),
            buffer: String::new(),
            header: String::new(),
            seq: String::new(),
        }
    }
}

impl<R: BufRead> Iterator for FastaIterator<R> {

    type Item = Sequence;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            self.buffer.clear();
            match self.reader.read_line(&mut self.buffer) {
                Ok(0) => {
                    if self.seq.len() > 0 {
                        let ret = Some(Sequence::new(&self.header, &self.seq));
                        self.seq.clear();
                        return ret;
                    }
                    return None;
                },
                Ok(_) => {
                    let line = self.buffer.trim();
                    if line.starts_with('>') {                  // --- It's a header!
                        if self.seq.len() > 0 {                      // --- we already have a sequence to return
                            let ret = Some(Sequence::new(&self.header, &self.seq));
                            self.header = self.buffer[1..].trim().to_owned();
                            self.seq.clear();
                            return ret;
                        } else {
                            self.header = self.buffer[1..].trim().to_owned();
                        }
                    } else {                                         // --- It's sequence
                        if line.len() > 0 { self.seq.push_str(line); }
                    }
                }
                Err(_) => return None,
            }
        }
    }
}


/// Iterator that reads sequences from a Stockholm-formatted buffer.
///
/// The Stockholm format is typically used to store a multiple alignment, e.g. used by HMMER, Pfam, and Rfam.
/// The first line of a Stockholm file (``.sto``, or ``.stk``) states the format and version identifier,
/// currently ``# STOCKHOLM 1.0``. The header is followed by mark-up lines beginning with ``#``.
/// These mark-up lines can annotate features of the alignment file (``#=GF``, generic per-file annotation),
/// or features of the aligned sequences (``#=GS``, generic per-sequence annotation).
/// The sequence alignment itself is a series of lines with sequence names (typically in the form name/start-end)
/// followed by a space and the aligned sequence. A line with two forward slashes (``//``) indicates the end of the alignment.
///
/// For the detailed description and a list of allowed annotations, see [Stockholm specification](https://sonnhammer.sbc.su.se/Stockholm.html)
///
pub struct StockholmIterator<R> {
    data: Vec<Sequence>,
    idx: usize,
    phantom: PhantomData<R>,
}

impl<R: BufRead> StockholmIterator<R> {
    pub fn new(stream: &mut R) -> Self {
        let data = StockholmIterator::from_stockholm_reader(stream);
        StockholmIterator { data: data, idx: 0, phantom: PhantomData, }
    }

    /// Read sequences in Stockholm format.
    ///
    /// Currently the function reads only the sequences; their annotations are not loaded
    ///
    /// For the detailed description of the format, see: `<https://sonnhammer.sbc.su.se/Stockholm.html>`
    pub fn from_stockholm_reader(reader: &mut R) -> Vec<Sequence> {

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
            // --- Remove header produced by the muscle program
            if buffer.starts_with("MUSCLE (") { continue }
            // --- first 8 spaces meand there is no seq-id => no sequence entry here, just markup
            if buffer.starts_with("        ") { continue }
            if ! buffer.starts_with('#') {      // --- the reader doesn't parse sequence annotations yet
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
}

impl<S> Iterator for StockholmIterator<S> where S: BufRead {

    type Item = Sequence;

    fn next(&mut self) -> Option<Self::Item> {
        if self.idx < self.data.len() {
            self.idx += 1;
            return Some(self.data[self.idx-1].clone());
        } else { return None; }
    }
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
/// use bioshell_seq::sequence::{Sequence, a3m_to_fasta};
/// use bioshell_seq::sequence::A3mConversionMode::RemoveSmallCaps;
/// // create two Sequence objects
/// let mut sequences = vec![Sequence::from_attrs(String::from("seq1"), "MTYKLILNGKTLKgeTTTEAVDAAT".as_bytes().to_vec()),
///         Sequence::from_attrs(String::from("seq2"), "HQYKLALNGKTLKaqTTTEAVDAAT".as_bytes().to_vec())];
///
/// a3m_to_fasta(&mut sequences, &RemoveSmallCaps);
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
                sequences[i] = Sequence::from_attrs(seq.description.clone(), seq_data);
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
/// use bioshell_seq::sequence::{Sequence, StockholmIterator, remove_gaps_by_sequence};
/// use std::io::{BufReader};
///
/// let alignment = "A0A6G1KYC8/4-30                 -------TK----QT---TW-EK--PA------
/// A0A6G1KYC8/54-83                -------TK----AT---AW-EM--PK-QM---
/// A0A2Z6MTB8/60-86                -------TR----QS---SW-EK--P-------
/// A0A2Z6MTB8/102-129              -------TQ----QS---TW-TI--PEE-----
/// H2XMA8/16-49                    -------TQ----RT---TW-QD--PR------
/// UPI000A31515F/122-163           -------TR----ES---AW-TK--PD------
/// UPI000A31515F/383-441           -------TL----ES---TW-EK--PQE-----
/// UPI00057AF938/507-540           -------TR----TT---TW-KH--PC------
/// UPI00138FB958/985-1021          -------TQ----QT---SW-LH--PVSQ----";
/// let mut reader = BufReader::new(alignment.as_bytes());
/// let mut sequences = StockholmIterator::from_stockholm_reader(&mut reader);
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
        let id = &sequences[j].description;
        let old_aa: &Vec<u8> = sequences[j].seq();
        let mut new_aa: Vec<u8> = vec![];
        for pi in &pos {
            new_aa.push(old_aa[*pi]);
        }
        sequences[j] = Sequence::from_attrs(id.clone(), new_aa);
    }
}

/// Counts residues of a given type in a [`Sequence`](Sequence)
///
/// # Examples
/// ```rust
/// use bioshell_seq::sequence::{Sequence, count_residue_type};
/// let seq_str = b"MTYKLILNGKTLKGETTTEAVDAATAEKVFKQY";
/// let sequence = Sequence::from_attrs(String::from("test_seq"), seq_str.to_vec());
/// assert_eq!(count_residue_type(&sequence,'T'),6);
/// ```
pub fn count_residue_type(sequence: &Sequence, res_type: char) -> usize {
    sequence.seq.iter().filter(|&c| *c == res_type as u8).count()
}

/// Counts identical residues between two sequences.
///
/// Matching gap symbols are not included in the count. Returns
/// [`AlignedSequencesOfDifferentLengths`](crate::SequenceError::AlignedSequencesOfDifferentLengths)
///  if the sequences differ by length.
///
/// # Example
/// ```rust
/// # use bioshell_seq::sequence::{Sequence, count_identical};
/// let si = Sequence::from_str("seq-1", "PERF");
/// let sj = Sequence::from_str("seq-2", "P-RV");
/// assert_eq!(count_identical(&si, &sj).unwrap(), 2);
/// ```
pub fn count_identical(si: &Sequence, sj: &Sequence) -> Result<usize, SequenceError> {
    if si.len() != sj.len() {
        return Err(SequenceError::AlignedSequencesOfDifferentLengths {
            length_expected: si.len(),
            length_found: sj.len(),
        });
    }

    return Ok(MSA::sum_identical(si, sj));
}


/// Length of a [`Sequence`](Sequence) excluding gaps
///
/// # Examples
/// ```rust
/// use bioshell_seq::sequence::{Sequence, len_ungapped};
/// let sequence = Sequence::from_str("test_seq", "P-RF");
/// assert_eq!(len_ungapped(&sequence),3);
/// let sequence = Sequence::from_str("test_seq", "__PERF_");
/// assert_eq!(len_ungapped(&sequence),4);
/// ```
pub fn len_ungapped(sequence: &Sequence) -> usize {
    sequence.seq.iter().filter(|&c| *c != b'-' && *c != b'_').count()
}

/// Removes gaps from a [`Sequence`](Sequence)
///
/// # Examples
/// ```rust
/// use bioshell_seq::sequence::{Sequence, remove_gaps};
/// let mut sequence = Sequence::from_str("test_seq", "P-RF_");
/// remove_gaps(&mut sequence);
/// assert_eq!(sequence.to_string(), String::from("PRF"));
/// ```
pub fn remove_gaps(sequence: &mut Sequence) {

    let new_seq: Vec<_>  = sequence.seq.iter()
        .filter(|&c| *c != b'-' && *c != b'_')
        .map(|&c| c)
        .collect();
    sequence.seq.clear();
    for c in new_seq.iter() { sequence.seq.push(*c)}
}


/// Makes an un-gapped copy of a given [`Sequence`](Sequence)
///
/// Returns a new [`Sequence`](Sequence) object that is an un-gapped copy of the argument. The new
/// [`Sequence`](Sequence) will inherit a description string from its source.
///
/// # Examples
/// ```rust
/// use bioshell_seq::sequence::{Sequence, clone_ungapped};
/// let mut sequence = Sequence::from_str("test_seq", "P-RF_");
/// let ungapped = clone_ungapped(&mut sequence);
/// assert_eq!(ungapped.to_string(), String::from("PRF"));
/// ```
pub fn clone_ungapped(sequence: &Sequence) -> Sequence {

    let mut s = sequence.clone();
    remove_gaps(&mut s);
    return s
}