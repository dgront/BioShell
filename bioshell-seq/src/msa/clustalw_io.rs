use std::collections::HashMap;
use std::io::BufRead;
use crate::msa::MSA;
use crate::sequence::Sequence;
use crate::SequenceError;

impl MSA {
    /// Reads a ClustalW/Clustal-formatted alignment and returns an `MSA`.
    ///
    /// The parser accepts the usual block-based Clustal format where the file starts
    /// with a header line containing the word `CLUSTAL` and the alignment is given
    /// in one or more blocks.
    /// Consensus lines (with `*`, `:` or `.`) are ignored. Sequence fragments are
    /// concatenated per sequence id preserving first-seen ordering.
    ///
    /// ```
    /// # use bioshell_seq::sequence::FastaIterator;
    /// # use bioshell_seq::SequenceError;
    /// # use std::io::BufReader;
    /// # fn main() -> Result<(), SequenceError> {
    /// use bioshell_seq::msa::MSA;
    /// let clustal_input = "
    /// CLUSTAL W (1.82) multiple sequence alignment
    ///
    /// seq1     A--TGCC
    /// seq2     AATTG-C
    ///          * ** *
    ///
    /// seq1     AGGTA
    /// seq2     AGG-A
    ///          *** *
    /// ";
    ///
    /// let mut reader = BufReader::new(clustal_input.as_bytes());
    /// let msa = MSA::from_clustalw_reader(&mut reader)?;
    /// assert_eq!(msa.n_seq(), 2);
    /// assert_eq!(msa.sequences()[0].id(), "seq1");
    /// assert_eq!(msa.sequences()[0].seq(), b"A--TGCCAGGTA");
    /// assert_eq!(msa.sequences()[1].id(), "seq2");
    /// assert_eq!(msa.sequences()[1].seq(), b"AATTG-CAGG-A");
    ///
    /// # Ok(())
    /// # }
    pub fn from_clustalw_reader<R: BufRead>(reader: &mut R) -> Result<Self, SequenceError> {
        let mut seq_fragments: HashMap<String, String> = HashMap::new();
        let mut seq_order: Vec<String> = Vec::new();
        let mut saw_header = false;

        for line in reader.lines() {
            let line = line?;
            let text = line.trim_end();

            if text.is_empty() { continue; }

            // Try to locate header (may be preceded by comments/blank lines)
            if !saw_header {
                // header line typically contains the token "CLUSTAL" (case-insensitive)
                if text.to_uppercase().contains("CLUSTAL") {
                    saw_header = true;
                    continue;
                } else {
                    // ignore leading comments starting with '#'
                    if text.starts_with('#') { continue; }
                    // otherwise keep searching for header
                    continue;
                }
            }

            // After header: skip separator/consensus lines that don't look like sequence rows.
            // A sequence row is expected to contain at least two whitespace-separated tokens: id and residue-block.
            let mut parts = text.split_whitespace();
            if let (Some(first), Some(second)) = (parts.next(), parts.next()) {
                // If the first token looks like a consensus token (e.g. "*" or ":"), skip the line
                if first.chars().all(|c| !c.is_alphanumeric() && c != '_' && c != '-') {
                    continue;
                }

                if !seq_fragments.contains_key(first) {
                    seq_order.push(first.to_string());
                }

                seq_fragments.entry(first.to_string())
                    .and_modify(|s| s.push_str(second))
                    .or_insert_with(|| second.to_string());
            }
        }

        if !saw_header {
            return Err(SequenceError::InvalidClustalWFormat{ line: "".to_string(), description: "missing CLUSTAL header".to_string() });
        }

        let mut sequences = Vec::with_capacity(seq_order.len());
        for seq_id in seq_order {
            let seq = seq_fragments.remove(&seq_id).unwrap_or_default();
            sequences.push(Sequence::new(&seq_id, &seq));
        }

        let msa = MSA::from_sequences(sequences)?;
        Ok(msa)
    }
}
