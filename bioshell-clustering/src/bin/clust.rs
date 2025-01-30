use std::env;
use clap::{Parser};
use log::{info};
use std::time::Instant;
use bioshell_clustering::errors::ClusteringError;
use bioshell_clustering::hierarchical::{balance_clustering_tree, retrieve_data, hierarchical_clustering, retrieve_clusters, retrieve_data_id, medoid_by_min_max, retrieve_outliers, DistanceMatrix};
use bioshell_clustering::hierarchical::strategies::{average_link, centroid_link, complete_link, median_link, single_link, wards_method};
use bioshell_io::{out_writer};

#[derive(Parser, Debug)]
#[clap(name = "clust")]
#[clap(about = "Performs hierarchical clustering analysis.", long_about = None)]
struct Args {
    #[clap(long, short='d')]
    /// input file with distances: TSV format
    infile: String,
    #[clap(long)]
    /// don't cluster the data, detect outliers instead; an outlier is an element for which no other elements is within a given distance cutoff
    #[clap(long)]
    detect_outliers: Option<f32>,
    /// use the single linkage clustering (the default strategy)
    #[clap(long, action)]
    single: bool,
    /// use the complete linkage clustering (rather than single_link)
    #[clap(long, action)]
    complete: bool,
    /// use the average linkage clustering (rather than single_link)
    #[clap(long, action)]
    average: bool,
    /// use the complete linkage clustering (rather than single_link)
    #[clap(long, action)]
    centroid: bool,
    /// use the centroid linkage clustering (rather than single_link)
    #[clap(long, action)]
    median: bool,
    /// use the median linkage  clustering (rather than single_link)
    #[clap(long, action)]
    wards: bool,
    /// writes clusters created by stopping the clustering at a given distance cutoff
    #[clap(short='c', long)]
    cutoff: Option<f32>,
    // /// writes the clustering tree in the Newick format
    // #[clap(long)]
    // newick: Option<String>,
    /// writes the distance matrix ordered by the clustering tree
    #[clap(long)]
    reorder_distances: Option<String>,
}

impl Args {
    pub(crate) fn clustering_strategy(&self) -> ClusteringStrategy {
        if self.complete {
            ClusteringStrategy::CompleteLink
        } else if self.average {
            ClusteringStrategy::AverageLink
        } else if self.wards {
            ClusteringStrategy::WardsMethod
        } else if self.median {
            ClusteringStrategy::AverageLink
        } else if self.centroid {
            ClusteringStrategy::WardsMethod
        } else {
            ClusteringStrategy::SingleLink
        }
    }
}

pub enum ClusteringStrategy {
    SingleLink,
    CompleteLink,
    AverageLink,
    MedianLink,
    CentroidLink,
    WardsMethod,
}


pub fn main() -> Result<(), ClusteringError> {

    unsafe {
        if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    }
    env_logger::init();
    let args = Args::parse();

    // ---------- read the matrix of sequence identity values from a TSV file ----------
    let distance_matrix = DistanceMatrix::from_tsv(&args.infile)?;
    let n_data = distance_matrix.n_elements();

    // ---------- detect outlier sequences; do not cluster -----------
    if let Some(cutoff) = args.detect_outliers {
        let distance_fn = |i: usize, j: usize| distance_matrix.distance(i, j);
        let outliers = retrieve_outliers(n_data, &distance_fn, cutoff);
        for outlier in outliers {
            println!("{}", distance_matrix.element_id(outlier));
        }
        return Ok(());
    }

    // ---------- cluster the sequences using the sequence identity matrix as distances ----------
    let start = Instant::now();
    let distance_fn = |i: usize, j: usize| distance_matrix.distance(i, j);
    let strategy = args.clustering_strategy();
    let mut clustering = match strategy {
        ClusteringStrategy::SingleLink => hierarchical_clustering(n_data, &distance_fn, &single_link),
        ClusteringStrategy::CompleteLink => hierarchical_clustering(n_data, &distance_fn, &complete_link),
        ClusteringStrategy::AverageLink => hierarchical_clustering(n_data, &distance_fn, &average_link),
        ClusteringStrategy::MedianLink => hierarchical_clustering(n_data, &distance_fn, &median_link),
        ClusteringStrategy::CentroidLink => hierarchical_clustering(n_data, &distance_fn, &centroid_link),
        ClusteringStrategy::WardsMethod => hierarchical_clustering(n_data, &distance_fn, &wards_method),
    };
    info!("{} sequences clustered in {:?}", n_data, start.elapsed());

    // ---------- retrieve actual clusters ----------
    if let Some(cutoff) = args.cutoff {
        let mut clusters = retrieve_clusters(&mut clustering, cutoff);
        clusters.sort_by(|a, b| a.value.cluster_size.cmp(&b.value.cluster_size));
        info!("{} clusters rertrieved for seq_id {:?}", clusters.len(), cutoff);
        for (i, cluster) in clusters.iter().enumerate() {
            let mut out_file = out_writer(&format!("cluster_{}-{}.dat", i, cluster.value.cluster_size), false);
            let medoid_idx = medoid_by_min_max(cluster, &distance_fn);
            writeln!(out_file, "# size: {}", cluster.value.cluster_size)?;
            writeln!(out_file, "# medoid: {}", distance_matrix.element_id(medoid_idx))?;
            let leaf_ids: Vec<usize> = retrieve_data_id(&cluster);
            for id in &leaf_ids {
                writeln!(out_file, "{}", distance_matrix.element_id(*id))?;
            }
            out_file.flush().unwrap();
        }
    }

    // ---------- retrieve the re-ordered distance matrix ----------
    if let Some(fname) = args.reorder_distances {
        balance_clustering_tree(&mut clustering, &distance_fn);
        let indexes: Vec<usize> =  (0..n_data).collect();
        let seq_order = retrieve_data(&clustering, &indexes);

        let mut out_file = out_writer(&fname, false);
        let mut k = 0;
        for i in &seq_order {
            let mut l = 0;
            for j in &seq_order {
                writeln!(out_file, "{}\t{}\t{:6.3}\t{}\t{}",
                         distance_matrix.element_id(*i), distance_matrix.element_id(*j),
                         distance_matrix.distance(*i, *j), k,l)?;
                l += 1;
            }
            writeln!(out_file)?;
            k += 1;
        }
    }

    Ok(())
}