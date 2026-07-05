use crate::sequence::parse_sequence_id;

/// Defines how a sequence label is created from a sequence's description string.
pub enum LabelStyle {
    /// returns the first recognized database identifier, the string will be truncated
    /// to the first ``n`` characters or padded with spaces to reach the length of ``n``
    FirstId { sort: bool, n: usize  },
    /// returns all the recognized database identifier, truncated to the first ``n`` characters
    FullId { sort: bool, n: usize },
    /// returns the original description string, truncated to the first ``n`` characters
    Description { n: usize },
}


/// Creates a string used to label a sequence based on its description line.
///
/// Such a label is intended to be a short but unique sequence identifier. It can be created in three different ways, depending on the ``style`` parameter:
/// 1. ``LabelStyle::FirstId { sort: bool }`` - returns the first recognized identifier in the description line.
/// If the ``sort`` flag is true, sorts them by priority as defined by [`SeqIdList::sort()`](crate::sequence::SeqIdList::sort) method.
/// 2. ``LabelStyle::FullId { sort: bool }`` - returns all recognized identifiers in the description line, separated by the '|' character.
/// If the ``sort`` flag is true, sorts them in the order defined by [`SeqIdList::sort()`](crate::sequence::SeqIdList::sort) method.
/// 3. ``Description`` - returns the first ``n`` characters of the description line; if ``n`` is 0 returns the entire description line.
///
/// # Examples
///```
/// use bioshell_seq::sequence::{Sequence, sequence_label};
/// use bioshell_seq::sequence::LabelStyle::{FirstId, FullId};
/// let desc = "sp|A0A009IHW8|ABTIR_ACIB9 2' cyclic ADP-D-ribose synthase [taxid=1310613]";
/// assert_eq!(sequence_label(desc, &FirstId {sort: false, n: 15}), "  sp|A0A009IHW8");
/// assert_eq!(sequence_label(desc, &FullId {sort: true, n:0}), "sp|A0A009IHW8|ABTIR_ACIB9|taxid=1310613");
/// ```
pub fn sequence_label(description: &str, style: &LabelStyle) -> String {
    let mut ids = parse_sequence_id(description);
    match style {
        LabelStyle::FirstId { sort, n } => {
            if *sort { ids.sort(); }
            let out = ids[0].to_string();
            if *n == 0 {
                out
            } else if out.len() > *n {
                out[0..*n].to_string()
            } else {
                format!("{:>width$}", out, width = *n)
            }
        }
        LabelStyle::FullId { sort, n } => {
            if *sort { ids.sort(); }
            let desc = ids.to_string();
            if *n == 0 {
                desc
            } else {
                let len = desc.len().min(*n);
                desc[0..len].to_string()
            }
        }
        LabelStyle::Description { n } => {
            let len = description.len().min(*n);
            description[0..len].to_string()
        }
    }
}
