use std::marker::PhantomData;
use crate::alignment::{AlignmentPath, AlignmentStep};
use crate::scoring::SimilarityScore;

/// Aligns two sequences using the Smith-Waterman local alignment algorithm with Gotoh matrices.
///
/// # Example
/// ```
/// use bioshell_seq::alignment::LocalAlignment;
/// use bioshell_seq::scoring::{SequenceSimilarityScore, SubstitutionMatrix, SubstitutionMatrixList};
/// let blosum62 = SubstitutionMatrix::load(SubstitutionMatrixList::BLOSUM62);
/// let scoring = SequenceSimilarityScore::for_strings("TSAILDSLGAEEIRAYLP", "MQRPILDSLGNPTAEEVKAFHW", blosum62);
/// let mut aligner = LocalAlignment::new(22);
/// aligner.align(&scoring, -10, -2);
/// assert_eq!(aligner.recent_score(), 40);
/// ```
#[allow(non_snake_case)]
pub struct LocalAlignment<T: SimilarityScore> {
    recent_score: i32,
    query_length: usize,
    tmplt_length: usize,

    best_i: usize,
    best_j: usize,
    best_state: u8, // 0 = H, 1 = E, 2 = F

    H: Vec<i32>,      // The alignment matrix - the row being calculated
    E: Vec<i32>,      // Gotoh's helper matrix for the alignments ending with a gap in the query - the row being calculated (horizontal move)
    F: Vec<i32>,      // Gotoh's helper matrix for the alignments ending with a gap in the template - the row being calculated (vertical move)
    H_prev: Vec<i32>, // The alignment matrix - the very previous row
    E_prev: Vec<i32>, // Gotoh's helper matrix for the alignments ending with a gap in the query - the very previous row
    F_prev: Vec<i32>, // Gotoh's helper matrix for the alignments ending with a gap in the template - the very previous row

    arrows: Vec<Vec<u8>>,  // traceback for H: 0=STOP, 1=E, 2=diag, 4=F
    e_trace: Vec<Vec<u8>>, // traceback for E: 0=STOP, 1=extend, 2=open
    f_trace: Vec<Vec<u8>>, // traceback for F: 0=STOP, 1=extend, 2=open

    phantom: PhantomData<T>,
}

impl<T: SimilarityScore> LocalAlignment<T> {

    /// Creates a new [`LocalAlignment`] instance with the specified maximum sequence length.
    #[allow(non_snake_case)]
    pub fn new(max_seq_length: usize) -> LocalAlignment<T> {
        let max_seq_length = max_seq_length + 1;
        LocalAlignment {
            recent_score: 0,
            query_length: 0,
            tmplt_length: 0,

            best_i: 0,
            best_j: 0,
            best_state: 0,

            H: vec![0; max_seq_length],
            E: vec![0; max_seq_length],
            F: vec![0; max_seq_length],
            H_prev: vec![0; max_seq_length],
            E_prev: vec![0; max_seq_length],
            F_prev: vec![0; max_seq_length],

            arrows: vec![vec![0u8; max_seq_length]; max_seq_length],
            e_trace: vec![vec![0u8; max_seq_length]; max_seq_length],
            f_trace: vec![vec![0u8; max_seq_length]; max_seq_length],

            phantom: PhantomData,
        }
    }

    /// Aligns two sequences using the Smith-Waterman local alignment algorithm with Gotoh matrices.
    #[allow(non_snake_case)]
    pub fn align(&mut self, scoring: &T, gap_open: i32, gap_extend: i32) -> i32 {

        let mut H = &mut self.H;
        let mut H_prev = &mut self.H_prev;
        let mut E = &mut self.E;
        let mut E_prev = &mut self.E_prev;
        let mut F = &mut self.F;
        let mut F_prev = &mut self.F_prev;

        self.query_length = scoring.query_length();
        self.tmplt_length = scoring.template_length();
        self.recent_score = 0;
        self.best_i = 0;
        self.best_j = 0;
        self.best_state = 0;

        for row in &mut self.arrows { row.fill(0u8); }
        for row in &mut self.e_trace { row.fill(0u8); }
        for row in &mut self.f_trace { row.fill(0u8); }

        // ---------- Initialize local alignment borders with zeros
        for j in 0..=self.tmplt_length {
            H[j] = 0;
            E[j] = 0;
            F[j] = 0;
        }

        for i in 1..=self.query_length {
            std::mem::swap(&mut H, &mut H_prev);
            std::mem::swap(&mut E, &mut E_prev);
            std::mem::swap(&mut F, &mut F_prev);

            H[0] = 0;
            E[0] = 0;
            F[0] = 0;

            for j in 1..=self.tmplt_length {
                // --- horizontal gap progression: E[i][j]
                let e_from_e = E[j - 1] + gap_extend;
                let e_from_h = H[j - 1] + gap_open;
                let e_from_f = F[j - 1] + gap_open;
                let e_open = e_from_h.max(e_from_f);

                if e_from_e > 0 && e_from_e >= e_open {
                    E[j] = e_from_e;
                    self.e_trace[i][j] = 1u8; // extend E
                } else if e_open > 0 {
                    E[j] = e_open;
                    self.e_trace[i][j] = 2u8; // open from H/F
                } else {
                    E[j] = 0;
                    self.e_trace[i][j] = 0u8; // STOP
                }

                // --- vertical gap progression: F[i][j]
                let f_from_f = F_prev[j] + gap_extend;
                let f_from_h = H_prev[j] + gap_open;
                let f_from_e = E_prev[j] + gap_open;
                let f_open = f_from_h.max(f_from_e);

                if f_from_f > 0 && f_from_f >= f_open {
                    F[j] = f_from_f;
                    self.f_trace[i][j] = 1u8; // extend F
                } else if f_open > 0 {
                    F[j] = f_open;
                    self.f_trace[i][j] = 2u8; // open from H/E
                } else {
                    F[j] = 0;
                    self.f_trace[i][j] = 0u8; // STOP
                }

                // --- matching move
                let h_diag = H_prev[j - 1] + scoring.score(i - 1, j - 1); // --- the score of a match at i,j

                let mut h_val = 0;
                let mut arrow_flag: u8 = 0; // 0 = STOP

                if h_diag > h_val {
                    h_val = h_diag;
                    arrow_flag = 2u8;
                } else if h_diag == h_val && h_val > 0 {
                    arrow_flag |= 2u8;
                }

                if E[j] > h_val {
                    h_val = E[j];
                    arrow_flag = 1u8;
                } else if E[j] == h_val && h_val > 0 {
                    arrow_flag |= 1u8;
                }

                if F[j] > h_val {
                    h_val = F[j];
                    arrow_flag = 4u8;
                } else if F[j] == h_val && h_val > 0 {
                    arrow_flag |= 4u8;
                }

                H[j] = h_val;
                self.arrows[i][j] = arrow_flag;

                if H[j] > self.recent_score {
                    self.recent_score = H[j];
                    self.best_i = i;
                    self.best_j = j;
                    self.best_state = 0;
                }
                if E[j] > self.recent_score {
                    self.recent_score = E[j];
                    self.best_i = i;
                    self.best_j = j;
                    self.best_state = 1;
                }
                if F[j] > self.recent_score {
                    self.recent_score = F[j];
                    self.best_i = i;
                    self.best_j = j;
                    self.best_state = 2;
                }
            }
        }

        self.recent_score
    }

    /// Backtraces the most recent alignment.
    ///
    ///Returns analignment path as well as the start point of the alignment in the query
    /// and template sequences, respectively.
    pub fn backtrace(&self) -> (AlignmentPath, usize, usize) {

        #[derive(Clone, Copy, Debug, PartialEq, Eq)]
        enum TraceState { H, E, F }

        let mut i = self.best_i;
        let mut j = self.best_j;
        let mut state = match self.best_state {
            0 => TraceState::H,
            1 => TraceState::E,
            2 => TraceState::F,
            _ => panic!("Internal error in local backtrace(): invalid state"),
        };

        let mut cigar: Vec<AlignmentStep> = vec![];

        loop {
            match state {
                TraceState::H => {
                    let a = self.arrows[i][j];
                    if a == 0u8 {
                        break;
                    }
                    if a & 2u8 != 0 {
                        cigar.push(AlignmentStep::Match);
                        i -= 1;
                        j -= 1;
                        state = TraceState::H;
                    } else if a & 1u8 != 0 {
                        state = TraceState::E;
                    } else if a & 4u8 != 0 {
                        state = TraceState::F;
                    } else {
                        panic!("Internal error in local backtrace(): invalid H traceback");
                    }
                }

                TraceState::E => {
                    let t = self.e_trace[i][j];
                    if t == 0u8 {
                        break;
                    }
                    cigar.push(AlignmentStep::Horizontal);
                    j -= 1;
                    state = if t == 1u8 { TraceState::E } else { TraceState::H };
                }

                TraceState::F => {
                    let t = self.f_trace[i][j];
                    if t == 0u8 {
                        break;
                    }
                    cigar.push(AlignmentStep::Vertical);
                    i -= 1;
                    state = if t == 1u8 { TraceState::F } else { TraceState::H };
                }
            }
        }

        cigar.reverse();
        (AlignmentPath::from_attrs(cigar), i, j)
    }

    /// Returns the score of the most recent alignment.
    pub fn recent_score(&self) -> i32 {
        self.recent_score
    }

    /// Returns the end point of the most recent alignment in the query and template sequences, respectively.
    pub fn recent_end_point(&self) -> (usize, usize) {
        (self.best_i, self.best_j)
    }
}