use std::io;
use std::io::Write;
use std::path::Path;
use bioshell_io::out_writer;
use crate::sequence::{first_word_of_description, Sequence};

/// A trait for types that can report on a `Sequence`, e.g. write it to a file.
pub trait SequenceReporter {
    /// Reports on the provided `Sequence`.
    fn report(&mut self, seq: &Sequence) -> io::Result<()>;
}

pub struct WriteFasta {
    line_width: usize,
    writer: Box<dyn Write>
}

impl WriteFasta {
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

pub struct SplitFasta {
    line_width: usize,
    path: Option<String>,
}

impl SplitFasta {
    pub fn new(path: Option<String>, line_width: usize) -> Self {
        Self { line_width, path }
    }
}

impl SequenceReporter for SplitFasta {

    fn report(&mut self, seq: &Sequence) -> io::Result<()> {
        let seq_id = first_word_of_description(seq.description());
        let file_name = format!("{seq_id}.fasta");
        let out_path_fname = match &self.path {
            None => file_name,
            Some(path) => Path::new(path).join(file_name).to_string_lossy().to_string(),
        };
        let mut out = out_writer(&out_path_fname, false);
        write!(out, "{:width$}\n", seq, width = self.line_width)
    }
}