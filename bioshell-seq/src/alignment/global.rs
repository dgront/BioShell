use std::marker::PhantomData;
use crate::alignment::{AlignmentPath, AlignmentStep};
use crate::scoring::{SimilarityScore};

/// Aligns two protein sequences using the Needleman-Wunsh algorithm with Gotoh matrices
///
/// ```
/// use bioshell_seq::alignment::GlobalAligner;
/// use bioshell_seq::scoring::{SequenceSimilarityScore, SubstitutionMatrix, SubstitutionMatrixList};
/// let blosum62 = SubstitutionMatrix::load(SubstitutionMatrixList::BLOSUM62);
/// let scoring = SequenceSimilarityScore::for_strings("AKLY", "AGKIY", blosum62);
/// let mut aligner = GlobalAligner::new(5);
/// aligner.align(&scoring, -10, -2);
/// println!("{}", aligner.recent_score());
/// ```
#[allow(non_snake_case)]
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
    arrows: Vec<Vec<u8>>,
    e_trace: Vec<Vec<u8>>,
    f_trace: Vec<Vec<u8>>,
    phantom: PhantomData<T>,
}

impl<T:SimilarityScore> GlobalAligner<T> {
    #[allow(non_snake_case)]
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
            arrows: vec![vec![0u8; max_seq_length]; max_seq_length],
            e_trace: vec![vec![0u8; max_seq_length]; max_seq_length],
            f_trace: vec![vec![0u8; max_seq_length]; max_seq_length],
            phantom: PhantomData
        }
    }

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

        for row in &mut self.e_trace { row.fill(0u8); }
        for row in &mut self.f_trace { row.fill(0u8); }

        let impossible_score = gap_open * self.max_length as i32;

        // ---------- Initialize for tail gaps penalized
        let mut tmp = gap_open;
        H[0] = 0;
        self.recent_score = tmp - gap_open; // = 0
        E[0] = 0;
        F[0] = 0;

        for j in 1..=self.tmplt_length {
            F[j] = impossible_score;
            H[j] = tmp;
            E[j] = tmp;
            self.arrows[0][j] = 1u8;   // H comes from E on the top border
            self.e_trace[0][j] = 1u8;  // E extends leftwards along the top border
            tmp += gap_extend;
        }

        let mut gap_started_in_q = gap_open;
        for i in 1..=self.query_length {
            std::mem::swap(&mut H, &mut H_prev);
            std::mem::swap(&mut E, &mut E_prev);
            std::mem::swap(&mut F, &mut F_prev);

            self.arrows[i][0] = 4u8;   // H comes from F on the left border
            self.f_trace[i][0] = 1u8;  // F extends upwards along the left border

            H[0] = gap_started_in_q;
            F[0] = gap_started_in_q;
            E[0] = impossible_score;

            for j in 1..=self.tmplt_length {
                // --- horizontal gap progression: E[i][j]
                let e_from_e = E[j - 1] + gap_extend;
                let e_from_h = H[j - 1] + gap_open;
                let e_from_f = F[j - 1] + gap_open;

                if e_from_e >= e_from_h && e_from_e >= e_from_f {
                    E[j] = e_from_e;
                    self.e_trace[i][j] = 1u8; // extend E
                } else {
                    E[j] = e_from_h.max(e_from_f);
                    self.e_trace[i][j] = 0u8; // open gap from H or F
                }

                // --- vertical gap progression: F[i][j]
                let f_from_f = F_prev[j] + gap_extend;
                let f_from_h = H_prev[j] + gap_open;
                let f_from_e = E_prev[j] + gap_open;

                if f_from_f >= f_from_h && f_from_f >= f_from_e {
                    F[j] = f_from_f;
                    self.f_trace[i][j] = 1u8; // extend F
                } else {
                    F[j] = f_from_h.max(f_from_e);
                    self.f_trace[i][j] = 0u8; // open gap from H or E
                }

                // --- matching move
                let h_diag = H_prev[j - 1] + scoring.score(i - 1, j - 1); // --- the score of a match at i,j
                H[j] = h_diag.max(E[j]).max(F[j]);

                let mut arrow_flag: u8 = 0;
                if H[j] == E[j] { arrow_flag += 1; }
                if H[j] == h_diag { arrow_flag += 2; }
                if H[j] == F[j] { arrow_flag += 4; }
                self.arrows[i][j] = arrow_flag;
            }
            gap_started_in_q += gap_extend;
        }

        self.recent_score = H[self.tmplt_length];
        H[self.tmplt_length]
    }
    pub fn backtrace(&self) -> AlignmentPath {

        #[derive(Clone, Copy, Debug, PartialEq, Eq)]
        enum TraceState { H, E, F }

        let mut i = self.query_length;
        let mut j = self.tmplt_length;
        let mut state = TraceState::H;
        let mut cigar: Vec<AlignmentStep> = vec![];

        while i > 0 || j > 0 {
            match state {
                TraceState::H => {
                    let a = self.arrows[i][j];

                    if a & 2 != 0 {
                        cigar.push(AlignmentStep::Match);
                        i -= 1;
                        j -= 1;
                        state = TraceState::H;
                    } else if a & 1 != 0 {
                        state = TraceState::E;
                    } else if a & 4 != 0 {
                        state = TraceState::F;
                    } else {
                        panic!("Internal error in backtrace(): invalid H traceback state");
                    }
                }

                TraceState::E => {
                    cigar.push(AlignmentStep::Horizontal);
                    j -= 1;

                    if self.e_trace[i][j + 1] == 1u8 {
                        state = TraceState::E; // continue horizontal gap
                    } else {
                        state = TraceState::H; // gap was opened from H/F, return to H
                    }
                }

                TraceState::F => {
                    cigar.push(AlignmentStep::Vertical);
                    i -= 1;

                    if self.f_trace[i + 1][j] == 1u8 {
                        state = TraceState::F; // continue vertical gap
                    } else {
                        state = TraceState::H; // gap was opened from H/E, return to H
                    }
                }
            }
        }

        cigar.reverse();
        AlignmentPath::from_attrs(cigar)
    }
    pub fn recent_score(&self) -> i32 { self.recent_score }
}
