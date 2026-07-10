use std::io::{BufRead, BufReader};
use std::fmt;
use log::{debug, info};
use crate::sequence::Sequence;
use crate::SequenceError;


/// Defines how sequences are parsed by the [`FastaIterator`] reader.
///
/// # Examples
///
/// ```rust
/// use bioshell_seq::sequence::FastaIterator;
/// # use bioshell_seq::SequenceError;
/// # fn main() -> Result<(), SequenceError> {
/// use bioshell_seq::sequence::FastaParsingMode::{Raw, CleanProtein};
/// let sequences: &str = "
/// > Heterobasidion annosum 172174
/// LNFCWVLGFSRYG (repeat) DQAELDSVVRTDRLPDFTDEAVLPYVTALMKKVLRWRPA";
/// let mut seqs = FastaIterator::new(sequences.as_bytes(), CleanProtein);
/// let s1 = seqs.next().ok_or(SequenceError::NoInputSequences)?;
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone, Copy, Default)]
pub enum FastaParsingMode {
    /// nucleic and protein sequences are assumed to be correct and loaded _as is_
    #[default]
    Raw,
    /// [`FastaIterator`] loads a protein sequence, possibly with gaps; other characters are removed
    CleanProtein,
    /// [`FastaIterator`] loads a protein sequence, gaps and the STOP codon marked as '*' are also allowed
    CleanProteinStop,
    /// [`FastaIterator`] loads a DNA/RNA sequence, possibly with gaps; other characters are removed
    CleanNucleic,
    /// user-provided function is called to every sequence line of the `fasta` input
    Custom(fn(&str, &mut String)),
}

impl fmt::Display for FastaParsingMode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let name = match self {
            Self::Raw => "raw",
            Self::CleanProtein => "clean-protein",
            FastaParsingMode::CleanProteinStop => "clean-protein-stop",
            FastaParsingMode::CleanNucleic => "clean-nucleic",
            FastaParsingMode::Custom(_) => "custom",
        };

        f.write_str(name)
    }
}


/// This enum is used internally to provide the [`FastaIterator`] with the actual parsing code
/// and at the same time to keep it private
enum FastaParserState {
    Raw,
    CleanProtein(AllowedCharsOnly),
    CleanProteinStop(AllowedCharsOnly),
    CleanNucleic(AllowedCharsOnly),
    Custom(fn(&str, &mut String)),
}


/// Iterator that provides sequences from a FASTA-formatted buffer.
///
/// This object iterates over a buffer without loading its whole content which allows processing
/// very large FASTA files.
///
/// # Examples
/// ```rust
/// use bioshell_seq::sequence::FastaIterator;
/// # use bioshell_seq::SequenceError;
/// # fn main() -> Result<(), SequenceError> {
/// use bioshell_seq::sequence::FastaParsingMode::Raw;
/// let sequences: &str = "> 1clf:A
/// AYKIADSCVSCGACASECPVNAISQGDSIFVIDADTCIDCGNCANVCPVGAPVQE
///
/// > 1dur:A
/// AYVINDSCIACGACKPECPVNCIQEGSIYAIDADSCIDCGSCASVCPVGAPNPED
///
/// > 1fca:A
/// AYVINEACISCGACEPECPVDAISQGGSRYVIDADTCIDCGACAGVCPVDAPVQA";
/// let seqs = FastaIterator::new(sequences.as_bytes(), Raw);
/// let seq_ids: Vec<String> = seqs.map(|r| r.map(|s| s.first_id(true).to_string())).collect::<Result<Vec<_>,_>>()?;
/// assert_eq!(seq_ids, vec!["1clf:A", "1dur:A", "1fca:A"]);
/// # Ok(())
/// # }
/// ```
pub struct FastaIterator<R> {
    reader: BufReader<R>,
    buffer: String,
    header: String,
    seq: String,
    mode: FastaParserState
}

impl<R: BufRead> FastaIterator<R> {
    pub fn new(stream: R, parsing_strategy: FastaParsingMode) -> Self {
        let mode = match parsing_strategy {
            FastaParsingMode::Raw => FastaParserState::Raw,
            FastaParsingMode::CleanProtein => {
                FastaParserState::CleanProtein(AllowedCharsOnly::new(PROTEIN_LETTERS))
            },
            FastaParsingMode::CleanProteinStop => {
                FastaParserState::CleanProteinStop(AllowedCharsOnly::new(PROTEIN_LETTERS_STOP))
            }
            FastaParsingMode::CleanNucleic => {
                FastaParserState::CleanNucleic(AllowedCharsOnly::new(NUCLEIC_RESIDUES))
            },
            FastaParsingMode::Custom(f) => FastaParserState::Custom(f),
        };

        info!("Parsing a .fasta stream with {parsing_strategy} method");

        FastaIterator {
            reader: BufReader::new(stream),
            buffer: String::new(),
            header: String::new(),
            seq: String::new(),
            mode
        }
    }
}

impl<R: BufRead> Iterator for FastaIterator<R> {

    type Item = Result<Sequence, SequenceError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            self.buffer.clear();
            match self.reader.read_line(&mut self.buffer) {
                Ok(0) => {
                    if self.seq.len() > 0 {
                        let ret = Sequence::new(&self.header, &self.seq);
                        self.seq.clear();
                        return Some(Ok(ret));
                    }
                    return None;
                },
                Ok(_) => {
                    let line = self.buffer.trim();
                    if line.starts_with('#') {
                        return Some(Err(SequenceError::InvalidFastaFormat{
                            line: line.to_string(),
                            description: "Fasta line must not start with '#' character".to_string()
                        }));
                    }
                    if line.starts_with('>') {                  // --- It's a header!
                        if self.seq.len() > 0 {                      // --- we already have a sequence to return
                            let ret = Sequence::new(&self.header, &self.seq);
                            self.header = self.buffer[1..].trim().to_owned();
                            self.seq.clear();
                            return Some(Ok(ret));
                        } else {
                            self.header = self.buffer[1..].trim().to_owned();
                        }
                    } else {                                         // --- It's sequence
                        if !line.is_empty() {
                            let before = self.seq.len();
                            let n_input_aa = line.len();
                            match &mut self.mode {
                                FastaParserState::Raw => {
                                    self.seq.push_str(line);
                                }
                                FastaParserState::CleanProtein(ref mut cleaner) => {
                                    cleaner.parse_line(line, &mut self.seq)
                                }
                                FastaParserState::CleanProteinStop(ref mut cleaner) => {
                                    cleaner.parse_line(line, &mut self.seq)
                                }
                                FastaParserState::CleanNucleic(ref mut cleaner) => {
                                    cleaner.parse_line(line, &mut self.seq)
                                }
                                FastaParserState::Custom(func) => { func(line, &mut self.seq) }
                            }
                            let n_accepted_aa = self.seq.len() - before;
                            if n_accepted_aa != n_input_aa {
                                let start = self.seq.len().saturating_sub(n_accepted_aa);
                                debug!("sequence fragment sanitized: received {n_input_aa}, stored {n_accepted_aa} bytes\nreceived: {line}\nrecorded: {}",&self.seq[start..]);
                            }
                        }
                    }
                }
                Err(err) => return Some(Err(SequenceError::Io(err))),
            }
        }
    }
}


#[derive(Debug, Clone, Default)]
struct AllowedCharsOnly {
    inside_parentheses: bool,
    allowed_letters: &'static [u8]
}

impl AllowedCharsOnly {
    pub fn new(allowed_letters: &'static [u8]) -> Self {
        Self {
            inside_parentheses: false,
            allowed_letters
        }
    }

    #[inline]
    pub fn parse_line(&mut self, line: &str, out: &mut String) {
        for b in line.bytes() {
            match b {
                b'(' => { self.inside_parentheses = true; }
                b')' => { self.inside_parentheses = false; }
                // Skip annotation content, e.g. "(gap)".
                _ if self.inside_parentheses => {}
                _ if self.allowed_letters.contains(&b) => {
                    out.push((b as char).to_ascii_uppercase());
                }
                // Skip spaces, digits, punctuation, coordinates, etc.
                _ => {}
            }
        }
    }
}

/// Defines letters allowed in an nucleic acid sequence, gaps are also allowed
const NUCLEIC_RESIDUES: &[u8] = b"cgmtu-_";

/// Defines letters allowed in an amino acid sequence, gaps are also allowed
static PROTEIN_LETTERS: &[u8] = b"ACDEFGHIKLMNPQRSTVWYBJOUXZ-_";

/// Defines letters allowed in an amino acid sequence, as well as '*' for the stop codon and gaps
static PROTEIN_LETTERS_STOP: &[u8] = b"ACDEFGHIKLMNPQRSTVWYBJOUXZ-_*";

