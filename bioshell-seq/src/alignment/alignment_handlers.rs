use std::fmt::Write;

use crate::alignment::{aligned_symbols, AlignmentPath};
use crate::sequence::Sequence;


/// Formats a pairwise alignment in the `pairwise` format.
///
/// Given the query and template sequences, an alignment path, and formatting parameters, this function returns a string representation of the alignment in a human-readable format.
/// `seq_name_width` is the width of the sequence name column, and `alignment_width` is the maximum number of alignment columns per block.
pub fn format_pairwise(query: &Sequence, template: &Sequence, path: &AlignmentPath, seq_name_width: usize, alignment_width: usize) -> String {

    assert!(alignment_width > 0, "alignment_width must be positive");

    let (aligned_q, aligned_t) = aligned_symbols(path, &query.as_u8(), &template.as_u8(), b'-');

    let q_name = query.description_n(seq_name_width);
    let t_name = template.description_n(seq_name_width);

    let mut q_pos = path.query_start();
    let mut t_pos = path.template_start();
    let mut out = String::new();

    const NUM_WIDTH: usize = 6;

    for (q, t) in aligned_q.chunks(alignment_width).zip(aligned_t.chunks(alignment_width)) {
        let q_count = q.iter().filter(|&&c| c != b'-').count();
        let t_count = t.iter().filter(|&&c| c != b'-').count();

        let q_from = q_pos + 1;
        let t_from = t_pos + 1;

        q_pos += q_count;
        t_pos += t_count;

        let middle: String = q.iter().zip(t.iter()).map(|(&a, &b)| { if a == b && a != b'-' { '|' } else { ' ' } }).collect();

        let q_str = std::str::from_utf8(q).unwrap();
        let t_str = std::str::from_utf8(t).unwrap();

        writeln!(out, "{q_name:<seq_name_width$} {q_from:>NUM_WIDTH$} {q_str} {q_pos}").unwrap();
        writeln!(out, "{:seq_name_width$} {:NUM_WIDTH$} {middle}", "", "").unwrap();
        writeln!(out, "{t_name:<seq_name_width$} {t_from:>NUM_WIDTH$} {t_str} {t_pos}").unwrap();
        writeln!(out).unwrap();
    }

    out
}
