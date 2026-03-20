use std::fmt;
use crate::sequence::Sequence;

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
        let line_width = f.width().unwrap_or(0);
        write!(f, "> {}\n{}\n", self.description(), self.to_string(line_width))
    }
}