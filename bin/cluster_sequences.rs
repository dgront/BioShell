use std::collections::HashMap;
use std::env;
use std::fmt::{Display, format};
use std::io::Write;
use clap::{Parser};
use log::{debug, error, info};
use std::time::Instant;
use rand::Rng;
use bioshell_clustering::hierarchical::{balance_clustering_tree, retrieve_data, hierarchical_clustering, retrieve_clusters, retrieve_data_id, medoid_by_min_max, retrieve_outliers};
use bioshell_clustering::hierarchical::strategies::{average_link, complete_link, single_link};
use bioshell_seq::sequence::{bucket_clustering, load_sequences, Sequence};
use bioshell_io::{can_create_file, out_writer};
use bioshell_seq::alignment::{align_all_pairs, AlignmentReporter, AlignmentStatistics};
use bioshell_seq::scoring::SubstitutionMatrixList;
use bioshell_seq::SequenceError;

#[derive(Parser, Debug)]
#[clap(name = "cluster_sequences", version, author)]
#[clap(about = "Cluster amino acid sequences by sequence identity", long_about = None)]
struct Args {
    /// input file in FASTA format
    infile: String,
    /// gap opening penalty
    #[clap(long, default_value = "-10", short='o')]
    open: i32,
    /// gap extension penalty
    #[clap(long, default_value = "-2", short='e')]
    extend: i32,
    /// don't cluster the sequences, detect outliers instead; an outlier is a sequence for which no other sequence is within a given sequence identity fraction
    #[clap(long)]
    detect_outliers: Option<f32>,
    /// use the single linkage clustering (the default strategy)
    #[clap(long, action)]
    single_link: bool,
    /// use the complete linkage clustering (rather than single_link)
    #[clap(long, action)]
    complete_link: bool,
    /// use the average linkage clustering (rather than single_link)
    #[clap(long, action)]
    average_link: bool,
    /// writes clusters created by stopping the clustering at a given sequence identity fraction
    #[clap(short='c', long)]
    identity_cutoff: Option<f32>,
    /// print the representative sequence for each cluster (i.e. the medoid)
    #[clap(short='m', long, action)]
    medoids: bool,
    /// clusters sequences into buckets at the given sequence identity fraction
    #[clap(short='b', long)]
    bucket_clustering: Option<f32>,
    /// prefix to add to the output files: cluster sequences and the medoids; use something like "unique_job_id_" to avoid overwriting previous results
    #[clap(long, default_value = "")]
    prefix: String,
    // /// writes the clustering tree in the Newick format
    // #[clap(long)]
    // newick: Option<String>,
    /// writes the input sequences reordered according to the clustering tree
    #[clap(long)]
    fasta: Option<String>,
    /// writes the distance matrix ordered by the clustering tree
    #[clap(long)]
    distance_matrix: Option<String>,
    /// length of a sequence name to print; longer names will be trimmed that size
    #[clap(long, short='w', default_value = "20")]
    name_width: usize,
    /// length of a sequence itself to print; use 0 to print the whole sequence in a single line
    #[clap(long, default_value = "80")]
    sequence_width: usize,
    /// be more verbose and log program actions on the screen
    #[clap(short, long, short='v')]
    verbose: bool
}


struct SequenceIdentityMatrix {
    n_sequences: usize,
    description_to_index: HashMap<String, usize>,
    similarity_matrix: Vec<Vec<f32>>,
}

impl SequenceIdentityMatrix {
    pub fn new(sequences: &[Sequence], name_width: usize) -> Result<SequenceIdentityMatrix, SequenceError> {
        let mut description_to_index: HashMap<String, usize> = HashMap::new();

        for (i, sequence) in sequences.iter().enumerate() {
            let outcome = description_to_index.insert(sequence.description().to_string(), i);
            if outcome.is_some() {
                return Err(SequenceError::IdenticalSequenceDescriptions {
                    description: sequence.description_n(name_width).to_string(),
                });
            }
        }

        Ok(SequenceIdentityMatrix {
            n_sequences: sequences.len(),
            description_to_index,
            similarity_matrix: vec![vec![0.0; sequences.len()]; sequences.len()],
        })
    }

    pub fn percent_identity(&self, i: usize, j:usize) -> f32 { self.similarity_matrix[i][j] }

    pub fn sequence_description(&self, i: usize) -> &str {
        self.description_to_index.iter().find(|(_desc, index)| *index == &i).unwrap().0
    }
}

impl Display for SequenceIdentityMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for i in 0..self.n_sequences {
            for j in 0..self.n_sequences {
                write!(f, "{} {} {:.2}\t", i, j, self.similarity_matrix[i][j])?;
            }
            writeln!(f)?;
        }
        Ok(())
    }
}

impl AlignmentReporter for SequenceIdentityMatrix {
    fn report(&mut self, aligned_query: &Sequence, aligned_template: &Sequence) {
        let stats = AlignmentStatistics::from_sequences(aligned_query, aligned_template, 32);
        let query_index = self.description_to_index[aligned_query.description()];
        let template_index = self.description_to_index[aligned_template.description()];
        self.similarity_matrix[query_index][template_index] = stats.percent_identity() as f32;
    }
}


pub fn main() -> Result<(), SequenceError> {

    let args = Args::parse();
    unsafe {
        if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
        if args.verbose { env::set_var("RUST_LOG", "debug"); }
    }
    env_logger::init();

    // ---------- load sequences ----------
    let sequences = load_sequences(&args.infile, "")?;

    if let Some(cutoff) = args.bucket_clustering {
        // ---------- create a random name of a dummy file to check if we can write
        let filename: String = rand::thread_rng()
            .sample_iter(rand::distributions::Alphanumeric)
            .take(10)
            .map(char::from)
            .collect();
        let prefix = &args.prefix;
        if !can_create_file(&format!("{}{}", prefix, filename)) {
            error!("Can't write with prefix {}", prefix);
            return Ok(());
        }
        let start = Instant::now();
        info!("Bucket clustering of {} sequences aligned at {}", sequences.len(), cutoff);
        let clusters = bucket_clustering(&sequences, cutoff);
        info!("{} sequences clustered in {:?}", sequences.len(), start.elapsed());
        for (i, cluster) in clusters.iter().enumerate() {
            let mut out_file = out_writer(&format!("{}cluster_{}-{}.fasta",
                                                   args.prefix, i, cluster.len()), false);
            for sequence in cluster.iter() {
                writeln!(out_file, "{:width$}", sequence, width = args.sequence_width)?;
            }
            out_file.flush().unwrap();
        }
        return Ok(());
    }

    // ---------- align all sequence pairs ----------
    let mut matrix_reporter = SequenceIdentityMatrix::new(&sequences, args.name_width)?;
    let start = Instant::now();
    align_all_pairs(&sequences, &sequences, SubstitutionMatrixList::BLOSUM62,
                    args.open, args.extend, true, &mut matrix_reporter);
    info!("{} sequences aligned in {:?}", sequences.len(), start.elapsed());

    // ---------- detect outlier sequences; do not cluster -----------
    if let Some(cutoff) = args.detect_outliers {
        let distance_fn = |i: usize, j: usize| 100.0 - matrix_reporter.percent_identity(i, j);
        let outliers = retrieve_outliers(sequences.len(), &distance_fn, 100.0 - cutoff);
        for outlier in outliers {
            println!("{}", sequences[outlier]);
        }
        return Ok(());
    }

    // ---------- cluster the sequences using the sequence identity matrix as distances ----------
    let start = Instant::now();
    let distance_fn = |i: usize, j: usize| 100.0 - matrix_reporter.percent_identity(i, j);
    let mut clustering = hierarchical_clustering(sequences.len(), distance_fn, &complete_link);
    info!("{} sequences clustered in {:?}", sequences.len(), start.elapsed());

    // ---------- retrieve actual clusters of sequences----------
    if let Some(cutoff) = args.identity_cutoff {
        let mut clusters = retrieve_clusters(&mut clustering, 100.0 - cutoff);
        clusters.sort_by(|a, b| a.value.cluster_size.cmp(&b.value.cluster_size));
        info!("{} clusters rertrieved for seq_id {:?}", clusters.len(), cutoff);
        for (i, cluster) in clusters.iter().enumerate() {
            let mut out_file = out_writer(&format!("{}cluster_{}-{}.fasta",
                    args.prefix, i, cluster.value.cluster_size), false);
            let leaf_ids: Vec<usize> = retrieve_data_id(&cluster);
            for id in &leaf_ids {
                let sequence = sequences.get(*id).unwrap();
                writeln!(out_file, "{:width$}", sequence, width = args.sequence_width)?;
            }
            out_file.flush().unwrap();
        }

        if args.medoids {
            for (i, cluster) in clusters.iter().enumerate() {
                let medoid_idx = medoid_by_min_max(cluster, &distance_fn);
                let mut out_file = out_writer(&format!("{}center_{}-{}.fasta",
                        args.prefix, i, cluster.value.cluster_size), false);
                writeln!(out_file, "{:width$}", &sequences[medoid_idx], width = args.sequence_width)?;
            }
        }
    }

    balance_clustering_tree(&mut clustering, &distance_fn);
    let indexes: Vec<usize> =  (0..sequences.len()).collect();
    let seq_order = retrieve_data(&clustering, &indexes);

    // ---------- retrieve the re-ordered distance matrix ----------
    if let Some(fname) = args.distance_matrix {
        let mut out_file = out_writer(&fname, false);
        let mut k = 0;
        for i in &seq_order{
            let mut l = 0;
            for j in &seq_order {
                writeln!(out_file, "{:4} {:4} {:6.3} {} {}", k, l, matrix_reporter.percent_identity(*i, *j),
                         matrix_reporter.sequence_description(*i), matrix_reporter.sequence_description(*j))?;
                l += 1;
            }
            writeln!(out_file)?;
            k += 1;
        }
    }

    // ---------- write the re-ordered sequences as FASTA ----------
    if let Some(fname) = args.fasta {
        let mut out_file = out_writer(&fname, false);
        for i in &seq_order {
            writeln!(out_file, "{}", sequences[*i])?;
        }
    }

    Ok(())
}