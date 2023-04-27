//! OPTICS clustering algorithm.
//!
use std::cmp::*;
use std::fmt;
use std::ops::{Index, Range};

use crate::{DistanceByIndex, CartesianPoints};
use crate::points_set::{Neighbor, NeighborsOf};


/// Defines the requirements for data points used in OPTICS clustering
///
/// Optics clustering accepts a points set that can compute distance between them and list its neighbors
pub trait OpticsPoints: DistanceByIndex + NeighborsOf {}

impl<T, D> OpticsPoints for CartesianPoints<T, D>
    where T: Index<usize, Output = f64>, D: Fn(&T, &T, usize) -> f64 {}

/// Provides the OPTICS (Ordering Points To Identify the Clustering Structure) clustering algorithm.
///
/// OPTICS is an algorithm for density based clustering. It looks for a  neighbourhood
/// of a certain radius `eps` must contain at least a minimal number of objects `min_points`.
/// OPTICS arranges input data in order, storing the core distance and a reasonable reachability
/// distance for each item.
pub struct Optics {
    /// The radius of a neighborhood.
    pub eps: f64,
    /// The minimum number of points required to form a dense region.
    pub min_points: usize,
    /// Data structure that holds the input points and can calculate the distance between them
    neighbors: Box<dyn OpticsPoints>,
    /// The total number of points to be clustered
    n: usize,
    /// The distance value used to mark unreachable points ("infinite" distance)
    undefined: f64,
    /// the output of the OPTICS algorithm: reachability distance for every point
    reachability: Vec<f64>,
    /// the output of the OPTICS algorithm: the order in which points were processed
    clustering_order: Vec<usize>,
    /// Flags saying whether a given point has been already processed or not
    processed: Vec<bool>,
    /// vector that serves as lazy priority queue of points to be processed
    seeds: Vec<usize>,
    /// Each cluster is defined as a range of the [`clustering_order`](Optics::clustering_order) vector
    clusters: Vec<Range<usize>>
}

impl Optics {

    /// Create clustering for a given set of parameters and an input data set.
    pub fn new(eps: f64, min_samples: usize, neighbors: Box<dyn OpticsPoints>) -> Self {
        let mut o = Self { eps, min_points: min_samples, neighbors, n:0, undefined: eps*10.0,
            reachability: vec![], clustering_order: vec![],
            processed: vec![], seeds: vec![], clusters: vec![] };
        o.run_clustering();

        return o;
    }

    /// Read-only access to the reacheability distance for each point; listed in the order defined
    /// by the [`clustering_order()`](Optics::clustering_order) method
    pub fn reacheability_distance(&self) -> &Vec<f64> { &self.reachability }

    /// Read-only access to the order in which points were clustered.
    /// The clustering order together with [`reacheability_distance()`](Optics::reacheability_distance) array
    /// provide so called *reacheability graph* of a OPTICS clustering process
    pub fn clustering_order(&self) -> &Vec<usize> { &self.clustering_order }

    /// Copy of clusters created by the clustering process.
    /// Each cluster is represented as a ``Vec<usize>`` vector if indexes pointing to the actuall data
    /// subjected to clustering
    pub fn clusters(&self) -> Vec<Vec<usize>> {
        let mut out: Vec<Vec<usize>> = Vec::new();
        for c_rg in &self.clusters {
            let cluster: Vec<usize> = c_rg.clone().map(|i| self.clustering_order[i]).collect();
            out.push(cluster);
        }
        return out;
    }

    /// Run the clustering
    fn run_clustering(&mut self) {
        self.n = self.neighbors.count_points();

        self.allocate();

        // --- loop over all points in the dataset
        for i in 0..self.n {
            if self.processed[i] { continue }      // --- skip already processed points
            // --- mark the point as processed, add it to the clustering order
            self.clustering_order.push(i);
            // --- find its neighbors
            let nb = self.neighbors.neighbors_of(i, self.eps);
            self.processed[i] = true;

            // --- skip the point if it's not a core point
            if nb.len() < self.min_points {continue}

            // --- ensure the list of seeds is empty (it should be actually)
            self.seeds.clear();
            // --- this will update reachability distances and insert core points to the list of seeds
            self.update_neighbors(&nb);
            while !self.seeds.is_empty() {
                // --- pop the point with currently lowest reachability distance value from the priority queue
                let pi = self.seeds[0];
                self.seeds[0] = self.seeds[self.seeds.len() - 1];
                self.seeds.truncate(self.seeds.len() - 1);
                self.seeds.sort_by(|lhs, rhs| self.reachability[*lhs].partial_cmp(&self.reachability[*rhs]).unwrap());
                // --- visit it
                self.clustering_order.push(pi);
                // --- find its neighbors
                let neighbrs = self.neighbors.neighbors_of(pi, self.eps);
                self.processed[pi] = true;
                if neighbrs.len() < self.min_points { continue; }     // --- the pi point is not a core point

                // --- update the neighbors
                self.update_neighbors(&neighbrs);
            }
        }

        assert_eq!(self.n, self.clustering_order.len());

        let mut start: usize = 0;
        for i in 1..self.n {
            let ci = self.clustering_order[i];
            if self.reachability[ci] == self.undefined {
                self.clusters.push(Range { start, end: i });
                start = i;
            }
        }
        self.clusters.push(Range { start, end: self.n });
    }

    fn update_neighbors(&mut self, neighbrs: &Vec<Neighbor>) {

        let core_dist = neighbrs[self.min_points - 1].d;
        for n in neighbrs {
            if self.processed[n.idx] {continue}
            let new_reach_dist = n.d.max(core_dist);
            if self.reachability[n.idx] == self.undefined { // n.idx was not inserted to seed yet
                self.seeds.push(n.idx);
            }
            if self.reachability[n.idx] > new_reach_dist {
                self.reachability[n.idx] = new_reach_dist; // set or update its reacheability distance
            }
        }
        self.seeds.sort_by(|lhs, rhs| self.reachability[*lhs].partial_cmp(&self.reachability[*rhs]).unwrap());
    }

    /// Allocate internal buffers for a given number of points subjected to clustering
    fn allocate(&mut self) {
        self.reachability = vec![self.undefined; self.n];   // epsilon * 10 is already the "infinite" distance that defines points that can't be reached
        self.processed = vec![false; self.n];
    }
}


// wrapped into a new type to be able to print it nicely
struct NeighborList(Vec<Neighbor>);

// prints the list of neighbors
impl fmt::Display for NeighborList {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.0.iter().fold(Ok(()), |result, pt| {
            result.and_then(|_| write!(f, "{}", pt))
        })
    }
}
