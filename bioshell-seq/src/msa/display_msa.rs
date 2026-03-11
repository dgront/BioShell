use std::fmt;
use crate::msa::MSA;

impl fmt::Display for MSA {
    /// Creates a `String` representation of an `MSA` - FASTA format
    /// # Examples
    ///
    /// Reads an `MSA` and writes in in `.fasta` format
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
        for seq in self.sequences() {
            write!(f, "> {}\n{}\n", seq.description(), seq.to_string(line_width))?
        }
        Ok(())
    }
}

