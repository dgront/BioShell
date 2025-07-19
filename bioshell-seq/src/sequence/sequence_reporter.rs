use std::io;
use std::io::Write;
use std::path::Path;
use bioshell_io::out_writer;
use crate::sequence::{parse_sequence_id, SeqIdList, Sequence};

/// A trait for types that can report on a `Sequence`, e.g. write it to a file.
///
pub trait SequenceReporter {
    /// Reports on the provided `Sequence`.
    fn report(&mut self, seq: &Sequence) -> io::Result<()>;
}

/// A [`SequenceReporter`](SequenceReporter) that appends a [`Sequence`](Sequence) to a file in FASTA format.
///
/// All the sequences reported by this struct are written to the same file.
pub struct WriteFasta {
    line_width: usize,
    writer: Box<dyn Write>
}

impl WriteFasta {
    /// Creates a new [`WriteFasta`](WriteFasta) object.
    ///
    /// If `file_name` is `None`, or `"stdout"` the output is written to stdout. Otherwise, a new file is created.
    /// If `if_append` is `true`, the output is appended to the file if it exists.
    /// Sequence is wrapped every `line_width` characters.
    pub fn new(file_name: Option<String>, line_width: usize, if_append: bool) -> Self {
        let out = match &file_name {
            None => out_writer("", true),
            Some(fname) =>  out_writer(fname, if_append)
        };
        Self { line_width, writer: out }
    }
}

impl SequenceReporter for WriteFasta {
    fn report(&mut self, seq: &Sequence) -> io::Result<()> {
        write!(self.writer, "{:width$}\n", seq, width = self.line_width)
    }
}

/// A [`SequenceReporter`](SequenceReporter) that writes each [`Sequence`](Sequence) to a *separate* file in FASTA format.
///
/// Each sequence is written to a separate file, named after the first word of its description.
pub struct SplitFasta {
    line_width: usize,
    path: Option<String>,
}

impl SplitFasta {
    /// Creates a new [`SplitFasta`](SplitFasta) object.
    ///
    /// The files are written to a directory specified by `path`. If `path` is `None`, the files are written to the current directory.
    /// Sequence is wrapped every `line_width` characters.
    pub fn new(path: Option<String>, line_width: usize) -> Self {
        Self { line_width, path }
    }
}

impl SequenceReporter for SplitFasta {

    fn report(&mut self, seq: &Sequence) -> io::Result<()> {
        let seq_id: SeqIdList = parse_sequence_id(seq.description());
        let file_name = format!("{}.fasta", seq_id.file_name());
        let out_path_fname = match &self.path {
            None => file_name,
            Some(path) => Path::new(path).join(file_name).to_string_lossy().to_string(),
        };
        let mut out = out_writer(&out_path_fname, false);
        write!(out, "{:width$}\n", seq, width = self.line_width)
    }
}