use std::marker::PhantomData;
use crate::alignment::{Aligner, AlignmentPath, AlignmentStep};
use crate::scoring::{SimilarityScore, SubstitutionMatrix};

macro_rules! max {
    ($x: expr) => ($x);
    ($x: expr, $($z: expr),+) => {{
        let y = max!($($z),*);
        if $x > y {
            $x
        } else {
            y
        }
    }}
}

/// Aligns two protein sequences using the Needleman-Wunsh algorithm with Gotoh matrices
///
/// ```
/// use bioshell_seq::alignment::GlobalAligner;
/// use bioshell_seq::scoring::{SequenceSimilarityScore, SubstitutionMatrix, SubstitutionMatrixList};
/// let blosum62 = SubstitutionMatrix::load(SubstitutionMatrixList::BLOSUM62);
/// let scoring = SequenceSimilarityScore::for_strings("AKLY", "AGKIY", blosum62);
/// let mut aligner = GlobalAligner::new(5);
/// aligner.align(scoring, -10, -2);
/// println!("{}", aligner.recent_score());
/// ```
pub struct GlobalAligner<T:SimilarityScore> {
    recent_score: i32,
    max_length: usize,
    query_length: usize,
    tmplt_length: usize,
    H: Vec<i32>, // The alignment matrix - the row being calculated
    E: Vec<i32>, // Gotoh's helper matrix for the alignments ending with a gap in the query - the row being calculated (horizontal move)
    F: Vec<i32>, // Gotoh's helper matrix for the alignments ending with a gap in the template - the row being calculated (vertical move)
    H_prev: Vec<i32>, // The alignment matrix - the very previous row
    E_prev: Vec<i32>, // Gotoh's helper matrix for the alignments ending with a gap in the query - the very previous row
    F_prev: Vec<i32>, // Gotoh's helper matrix for the alignments ending with a gap in the template - the very previous row
    arrows: Vec<Vec<u16>>,
    E_arrows: Vec<Vec<u16>>,
    F_arrows: Vec<Vec<u16>>,
    phantom: PhantomData<T>,
}

impl<T:SimilarityScore> GlobalAligner<T> {
    pub fn new(max_seq_length: usize) -> GlobalAligner<T> {
        let max_seq_length = max_seq_length + 1;
        GlobalAligner {
            recent_score: 0,
            max_length: max_seq_length,
            query_length: 0,
            tmplt_length: 0,
            H: vec![0; max_seq_length],
            E: vec![0; max_seq_length],
            F: vec![0; max_seq_length],
            H_prev: vec![0; max_seq_length],
            E_prev: vec![0; max_seq_length],
            F_prev: vec![0; max_seq_length],
            arrows: vec![vec![0u16; max_seq_length]; max_seq_length],
            E_arrows: vec![vec![0u16; max_seq_length]; max_seq_length],
            F_arrows: vec![vec![0u16; max_seq_length]; max_seq_length],
            phantom: PhantomData
        }
    }

    pub fn align(&mut self, scoring: &T, gap_open: i32, gap_extend: i32) -> i32 {

        self.query_length = scoring.query_length();
        self.tmplt_length = scoring.template_length();
        for e in &mut self.E_arrows { e.fill(1u16); }    // open a new gap everywhere
        for f in &mut self.F_arrows { f.fill(1u16); }    // open a new gap everywhere
        let impossible_score = gap_open * self.max_length as i32;

        // ---------- Initialize for tail gaps penalized
        let mut tmp = gap_open;
        self.H[0] = 0;
        self.recent_score = 0;
        self.E[0] = 0;
        self.F[0] = 0;
        for j in 1..self.tmplt_length + 1 {
            self.F[j] = impossible_score;
            self.H[j] = tmp;
            self.E[j] = tmp;
            self.E_arrows[0][j] = j as u16;
            self.arrows[0][j] = 1;   // horizontal
            tmp += gap_extend;
        }

        let mut gap_started_in_q = gap_open;
        for i in 1..self.query_length + 1 {
            std::mem::swap(&mut self.H, &mut self.H_prev);
            std::mem::swap(&mut self.E, &mut self.E_prev);
            std::mem::swap(&mut self.F, &mut self.F_prev);
            self.arrows[i][0] = 4;   // vertical gap
            self.F_arrows[i][0] = i as u16;   // horizontal
            self.H[0] = gap_started_in_q;
            self.F[0] = gap_started_in_q;
            self.E[0] = impossible_score;
            for j in 1..self.tmplt_length + 1 {
                // --- horizontal gap progression
                self.E[j] = max!(self.E[j - 1] + gap_extend, self.H[j - 1] + gap_open, self.F[j - 1] + gap_open);
                if self.E[j] - self.E[j - 1] == gap_extend { self.E_arrows[i][j] = self.E_arrows[i][j-1] + 1; }
                // --- vertical gap progression
                self.F[j] = max!(self.F_prev[j] + gap_extend, self.H_prev[j] + gap_open, self.E_prev[j] + gap_open);
                if self.F[j] - self.F_prev[j] == gap_extend { self.F_arrows[i][j] = self.F_arrows[i - 1][j] + 1; }
                // --- matching move
                let h = self.H_prev[j - 1] + scoring.score(i - 1, j - 1); // --- the score of a match at i,j
                self.H[j] = max!(h, self.E[j], self.F[j]);
                let mut arrow_flag: u8 = 0;
                if self.H[j] == self.E[j] { arrow_flag += 1; }
                if self.H[j] == h { arrow_flag += 2; }
                if self.H[j] == self.F[j] { arrow_flag += 4; }
                self.arrows[i][j] = arrow_flag as u16;
            }
            gap_started_in_q += gap_extend;
        }
        self.recent_score = self.H[self.tmplt_length];

        return self.H[self.tmplt_length];
    }

    pub fn backtrace(&self) -> AlignmentPath {

        let mut i = self.query_length;
        let mut j = self.tmplt_length;
        let mut cigar: Vec<AlignmentStep> = vec![];
        while i > 0 || j > 0 {
            let a = self.arrows[i][j];
            if a & 2 != 0 {
                cigar.push(AlignmentStep::Match);
                i -= 1;
                j -= 1;
                continue;
            }
            if a & 1 != 0 {
                let n_gaps = self.E_arrows[i][j];
                for kk in 0..n_gaps {
                    cigar.push(AlignmentStep::Horizontal);
                    j -= 1;
                }
                continue;
            }
            if a & 4 != 0 {
                let n_gaps = self.F_arrows[i][j];
                for kk in 0..n_gaps {
                    cigar.push(AlignmentStep::Vertical);
                    i -= 1;
                }
                continue;
            }
            panic!("Internal error in backtrace()");
        }
        cigar.reverse();

        return AlignmentPath::from_attrs(cigar);
    }

    pub fn recent_score(&self) -> i32 { self.recent_score }
}
