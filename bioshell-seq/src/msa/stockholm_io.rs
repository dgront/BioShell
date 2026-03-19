use std::collections::HashMap;
use std::fmt;
use std::io::BufRead;
use crate::msa::{GcEntryType, GfEntryType, GrEntryType, GsEntryType, MSA, StockholmMSA};
use crate::sequence::Sequence;
use crate::SequenceError;

const MAX_DESCRIPTION_LENGTH: usize = 40;

impl fmt::Display for StockholmMSA {
    /// Creates a `String` representation of an `MSA` - Stockholm format
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {

        write!(f, "# STOCKHOLM 1.0\n")?;
        let block_length = f.width().unwrap_or(0);

        let (block_length, n_blocks) =
            if block_length == 0 { (self.len(), 1) } else { (block_length, self.len() / block_length + 1) };

        let descriptions: Vec<_> = self.sequences().iter().map(|s| s.description()).collect();
        let mut dsc_length = descriptions.iter().map(|dsc| dsc.len()).max().unwrap_or(0);
        dsc_length = dsc_length.min(MAX_DESCRIPTION_LENGTH);
        for i_block in 0..n_blocks {
            for i_seq in 0..descriptions.len() {
                let i_seq_str =self.sequences()[i_seq].seq();
                let end = ((i_block + 1) * block_length).min(i_seq_str.len());
                let seq_part = String::from_utf8_lossy(&i_seq_str[i_block* block_length..end]);
                write!(f, "{:<dsc_length$}\t{}\n", descriptions[i_seq], seq_part)?;
            }
        }
        Ok(())
    }
}

impl StockholmMSA {
    /// Reads a Stockholm alignment, processing also annotations such as `#=GF` or `#=GS`
    ///
    /// # Example
    /// ```
    /// # use std::fs::File;
    /// # use std::io::BufReader;
    /// # use bioshell_seq::msa::StockholmMSA;
    /// # use bioshell_seq::SequenceError;
    /// # fn main() -> Result<(), SequenceError> {
    /// let file = File::open("tests/test_files/4Fe-4S-example.sto")?;
    /// let mut sto_reader = BufReader::new(file);
    /// let msa = StockholmMSA::from_stockholm_reader(&mut sto_reader)?;
    /// assert_eq!(msa.n_seq(), 6);
    /// Ok(())
    /// # }
    /// ```
    pub fn from_stockholm_reader<R: BufRead>(reader: &mut R) -> Result<Self, SequenceError> {
        let mut gf: HashMap<GfEntryType, Vec<String>> = HashMap::new();
        let mut gs: HashMap<String, HashMap<GsEntryType, Vec<String>>> = HashMap::new();
        let mut gc: HashMap<GcEntryType, String> = HashMap::new();
        let mut gr: HashMap<String, HashMap<GrEntryType, String>> = HashMap::new();

        // Sequence text accumulated per sequence id.
        let mut seq_fragments: HashMap<String, String> = HashMap::new();

        // Preserve first-seen order of sequences.
        let mut seq_order: Vec<String> = Vec::new();

        // Prefer Stockholm GS DE as description if present.
        let mut descriptions: HashMap<String, String> = HashMap::new();

        let mut saw_header = false;

        for line in reader.lines() {
            let line = line?;
            let text = line.trim();

            if text.is_empty() { continue; }

            if text == "# STOCKHOLM 1.0" {
                saw_header = true;
                continue;
            }

            if text == "//" { break; }

            if let Some(rest) = text.strip_prefix("#=GF ") {
                let mut parts = rest.splitn(2, char::is_whitespace);
                if let (Some(feature), Some(value)) = (parts.next(), parts.next()) {
                    gf.entry(GfEntryType::from(feature))
                        .or_default()
                        .push(value.trim().to_string());
                }
                continue;
            }

            if let Some(rest) = text.strip_prefix("#=GS ") {
                let mut parts = rest.splitn(3, char::is_whitespace);
                if let (Some(seq_id), Some(feature), Some(value)) =
                    (parts.next(), parts.next(), parts.next()) {
                    let value = value.trim().to_string();
                    let key = GsEntryType::from(feature);

                    gs.entry(seq_id.to_string())
                        .or_default()
                        .entry(key.clone())
                        .or_default()
                        .push(value.clone());

                    if key == GsEntryType::DE {
                        descriptions
                            .entry(seq_id.to_string())
                            .and_modify(|s| {
                                s.push(' ');
                                s.push_str(&value);
                            })
                            .or_insert(value);
                    }
                }
                continue;
            }

            if let Some(rest) = text.strip_prefix("#=GC ") {
                let mut parts = rest.splitn(2, char::is_whitespace);
                if let (Some(feature), Some(value)) = (parts.next(), parts.next()) {
                    gc.entry(GcEntryType::from(feature))
                        .and_modify(|s| s.push_str(value.trim()))
                        .or_insert_with(|| value.trim().to_string());
                }
                continue;
            }

            if let Some(rest) = text.strip_prefix("#=GR ") {
                let mut parts = rest.splitn(3, char::is_whitespace);
                if let (Some(seq_id), Some(feature), Some(value)) =
                    (parts.next(), parts.next(), parts.next())
                {
                    gr.entry(seq_id.to_string())
                        .or_default()
                        .entry(GrEntryType::from(feature))
                        .and_modify(|s| s.push_str(value.trim()))
                        .or_insert_with(|| value.trim().to_string());
                }
                continue;
            }


            // Regular Stockholm sequence row:
            // <seqname><whitespace><aligned-sequence>
            let mut parts = text.split_whitespace();
            if let (Some(seq_id), Some(seq_part)) = (parts.next(), parts.next()) {
                if !seq_fragments.contains_key(seq_id) {
                    seq_order.push(seq_id.to_string());
                }
                seq_fragments
                    .entry(seq_id.to_string())
                    .and_modify(|s| s.push_str(seq_part))
                    .or_insert_with(|| seq_part.to_string());
            }
        }

        if !saw_header {
            return Err(SequenceError::InvalidStockholmFormat{ line: "".to_string(), description: "missing Stockholm header '# STOCKHOLM 1.0'".to_string() });
        }

        let mut sequences = Vec::with_capacity(seq_order.len());
        for seq_id in seq_order {
            let seq = seq_fragments.remove(&seq_id).unwrap_or_default();
            // let desc = descriptions.remove(&seq_id).unwrap_or_default();

            // Adjust this line if your Sequence constructor differs.
            sequences.push(Sequence::new(&seq_id, &seq));
        }

        let msa = MSA::from_sequences(sequences)?;

        Ok(Self { msa, gf, gs, gc, gr })
    }
}
