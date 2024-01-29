use std::fmt::{Display, Formatter};
use std::fs::File;
use std::io::{BufRead, BufReader};
use crate::errors::ScoringError;
use crate::errors::ScoringError::{CantParseNCBIEntry, FileNotFound, IncorrectNCBIFormat, ReadingError};

/// Lists substitution matrices shipped with [`bioshell-seq`](crate) crate.
///
/// These matrices can be loaded as shown below:
/// ```
/// use bioshell_seq::scoring::{SubstitutionMatrix, SubstitutionMatrixList};
/// let blosum80 = SubstitutionMatrix::load(SubstitutionMatrixList::BLOSUM80);
/// assert_eq!(blosum80.score_by_aa('C' as u8, 'C' as u8), 13);
/// ```
pub enum SubstitutionMatrixList {
    BLOSUM45,
    BLOSUM80,
    PAM250,
    PAM70,
    BLOSUM62,
    PAM120,
    PAM30,
}


const N_COLUMNS: usize = 21;

/// Holds an amino acid substitution matrix (aka similarity matrix)
///
/// Such a matrix can be loaded from an external file in the NCBI format. See the [crate::scoring] documentation for an example.
pub struct SubstitutionMatrix {
    pub(crate) score: [[i32; N_COLUMNS]; N_COLUMNS],
    pub(crate) aa_indexes:[u8; 255],
    pub(crate) aa_codes:[u8; 21]
}

impl SubstitutionMatrix {
    pub(crate) fn new() -> SubstitutionMatrix {
        SubstitutionMatrix{ score: [[0; N_COLUMNS]; N_COLUMNS], aa_indexes: [0; 255], aa_codes: ['*' as u8; 21] }
    }

    /// Loads a [SubstitutionMatrix] that has been defined within the BioShell package
    ///
    /// This [`bioshell-seq`](bioshell-seq) crate provides a staple collection of substitution matrices,
    /// listed in the  [`SubstitutionMatrixList`](SubstitutionMatrixList) enum.
    ///
    /// # Example
    /// ```
    /// use bioshell_seq::scoring::{SubstitutionMatrix, SubstitutionMatrixList};
    /// let blosum80 = SubstitutionMatrix::load(SubstitutionMatrixList::BLOSUM80);
    /// assert_eq!(blosum80.score_by_aa('C' as u8, 'C' as u8), 13);
    /// ```
    pub fn load(matrix_name: SubstitutionMatrixList) -> SubstitutionMatrix {
        let data = match matrix_name {
            SubstitutionMatrixList::BLOSUM45 => { include_str!("../../data/BLOSUM45") }
            SubstitutionMatrixList::BLOSUM80 => { include_str!("../../data/BLOSUM80") }
            SubstitutionMatrixList::PAM250 => { include_str!("../../data/PAM250") }
            SubstitutionMatrixList::PAM70 => { include_str!("../../data/PAM70") }
            SubstitutionMatrixList::BLOSUM62 => { include_str!("../../data/BLOSUM62") }
            SubstitutionMatrixList::PAM120 => { include_str!("../../data/PAM120") }
            SubstitutionMatrixList::PAM30 => { include_str!("../../data/PAM30") }
        };
        return SubstitutionMatrix::ncbi_matrix_from_buffer(BufReader::new(data.as_bytes())).unwrap();
    }

    #[inline(always)]
    /// Returns index for an amino acid type
    ///
    /// The index points to a column (or a row) of a substitution matrix that holds values for that
    /// amino acid; E.g. ``aa_index('A' as u8)`` should return 0.
    pub fn aa_index(&self, aa_letter: u8) -> u8 { self.aa_indexes[aa_letter as usize] }

    #[inline(always)]
    /// Provides the score for a given pair of amino acid indexes according to this [SubstitutionMatrix].
    ///
    /// Amino acid indexes may be obtained by a [aa_index()](SubstitutionMatrix::aa_index()) method call
    pub fn score_by_index(&self, aa_index_i: u8, aa_index_j: u8) -> i32 {
        self.score[aa_index_i as usize][aa_index_j as usize]
    }

    #[inline(always)]
    /// Provides the score for a given pair of amino acids according to this [SubstitutionMatrix].
    ///
    /// Amino acid must be specified by their single-letter codes, such as ``'W'`` or ``'Q'``
    pub fn score_by_aa(&self, aa_letter_i: u8, aa_letter_j: u8) -> i32 {
        self.score_by_index(self.aa_index(aa_letter_i), self.aa_index(aa_letter_j))
    }

    /// Loads a [SubstitutionMatrix] from data in the NCBI format.
    ///
    /// Entries corresponding to ``'B'``, ``'J'``, ``'Z'`` and ``'*'`` symbols are not loaded.
    pub fn ncbi_matrix_from_buffer<R: BufRead>(reader: R) -> Result<SubstitutionMatrix, ScoringError>  {

        let mut m = SubstitutionMatrix::new();
        let mut i = 0_usize;
        for line in reader.lines() {
            let line = match line {
                Ok(l) => { l }
                Err(_) => {return Err(ReadingError)}
            };
            if line.starts_with("#") || line.starts_with(" ") { continue; }
            let values: Vec<&str> = line.split_whitespace().collect();
            let n_values = values.len();
            if n_values < 23 { return Err(IncorrectNCBIFormat{ line: line.clone() }) }
            let char_i = values[0].as_bytes()[0];
            m.aa_indexes[char_i as usize] = i as u8;
            m.aa_codes[i] = char_i;
            for j in 1..21 {
                m.score[i][j - 1] = match values[j].parse::<i32>() {
                    Ok(val) => { val }
                    Err(_) => { return Err(CantParseNCBIEntry{ line: line.clone(), value: values[j].to_string() }) }
                };
                m.score[j - 1][i] = m.score[i][j - 1]
            }
            // n_values - 2 is the index of 'X' column
            m.score[i][20] = match values[n_values - 2].parse::<i32>() {
                Ok(val) => { val }
                Err(_) => { return Err(CantParseNCBIEntry{ line: line.clone(), value: values[n_values - 2].to_string() }) }
            };
            m.score[20][i] = m.score[i][20];
            i += 1;
            if i==20 { break }
        }
        m.aa_indexes['X' as usize] = 20u8;
        m.aa_codes[20] = 'X' as u8;
        m.score[20][20] = -1;

        return Ok(m);
    }

    /// Loads a [SubstitutionMatrix] from a file in the NCBI format.
    ///
    /// This method simply opens the ``file_name`` file for reading and calls
    /// [ncbi_matrix_from_buffer()](SubstitutionMatrix::ncbi_matrix_from_buffer()).
    /// An example of a NCBI-formatted matrix may be found [here](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/BLOSUM62)
    pub fn ncbi_matrix_from_file(file_name: &str) -> Result<SubstitutionMatrix, ScoringError> {

        let file = match File::open(file_name) {
            Ok(f) => { f }
            Err(_) => { return Err(FileNotFound{ file_name: file_name.to_string() })}
        };
        let reader = BufReader::new(file);
        return SubstitutionMatrix::ncbi_matrix_from_buffer(reader);
    }
}

impl Display for SubstitutionMatrix {
    /// Displays a [SubstitutionMatrix] as a nice table
    ///
    /// **Note**, that the printed output only resembles the NCBI format! The [SubstitutionMatrix] struct
    /// doesn't use similarity scores for ``'B'``, ``'J'``, ``'Z'`` and ``'*'`` symbols. These
    /// values are not loaded from a respective NCBI file and are not shown in the output produced by this method
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f ,"# Entries for a matrix in the BioShell internal format\n")?;
        write!(f ,"   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  X\n")?;

        for i in 0..21 {
            write!(f, "{}", self.aa_codes[i] as char)?;
            for j in 0..21 {
                let score = self.score_by_aa(self.aa_codes[i], self.aa_codes[j]);
                write!(f, "{:3}", score)?;
            }
            write!(f,"\n")?;
        }
        Ok(())
    }
}
