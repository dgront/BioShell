use std::cmp::*;
use std::fmt;
use std::ops::{Range};

#[derive(Clone)]
pub struct PointAndDistance {
    pub idx: usize,
    pub d:f64
}

impl fmt::Display for PointAndDistance {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({}:{:4.3}) ", self.idx, self.d)
    }
}

pub trait Neighbors {

    fn distance(&self, i:usize, j:usize) -> f64;

    fn neighbors_of(&self, idx: usize, cutoff: f64) -> Vec<PointAndDistance>;

    fn count_points(&self) -> usize;
}

pub struct EuclideanNeighbors {
    datapoints: Vec<Vec<f64>>
}

impl EuclideanNeighbors {
    pub fn new(data: Vec<Vec<f64>>) -> EuclideanNeighbors {
        EuclideanNeighbors{datapoints: data}
    }
}


impl Neighbors for EuclideanNeighbors {

    fn count_points(&self) -> usize { self.datapoints.len() }

    fn distance(&self, i:usize, j:usize) -> f64 {
        let mut d: f64 = 0.0;
        let pi: &Vec<f64> = &self.datapoints[i];
        let pj: &Vec<f64> = &self.datapoints[j];
        for i in 0..pi.len() {
            let t = pi[i] - pj[i];
            d += t * t;
        }
        d.sqrt()
    }

    /// Finds neighbors of a given point.
    /// Collects all the points that are within a cutoff to the ``idx`` point; returns the results
    /// as a vector of [`PointAndDistance`](Optics::PointAndDistance) structs where the distance values are also stored
    fn neighbors_of(&self, idx:usize, cutoff: f64) -> Vec<PointAndDistance> {
        let mut out: Vec<PointAndDistance> = vec![];
        for i in 0..self.datapoints.len() {
            let d = self.distance(idx, i);
            if d <= cutoff {
                out.push(PointAndDistance { idx: i, d });
            }
        }
        // --- sort points by distance
        out.sort_by(|k, l| k.d.partial_cmp(&l.d).unwrap());
        return out;
    }
}

/// Provides the OPTICS (ordering points to identify the clustering structure) clustering algorithm.
pub struct Optics {
    /// The radius of a neighborhood.
    pub eps: f64,
    /// The minimum number of points required to form a dense region.
    pub min_points: usize,
    /// Data structure that holds the input points and can calculate the distance between them
    neighbors: Box<dyn Neighbors>,
    /// The total number of points to be clustered
    n: usize,
    /// The distance value used to mark unreachable points ("infinite" distance)
    undefined: f64,
    /// the output of the OPTICS algorithm: reacheability distance for every point
    reacheability: Vec<f64>,
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

    pub fn new(eps: f64, min_samples: usize, neighbors: Box<dyn Neighbors>) -> Self {
        let mut o = Self { eps, min_points: min_samples, neighbors, n:0, undefined: eps*10.0,
            reacheability: vec![], clustering_order: vec![],
            processed: vec![], seeds: vec![], clusters: vec![] };
        o.run_clustering();

        return o;
    }

    /// Read-only access to the reacheability distance for each point; listed in the order defined
    /// by the [`clustering_order()`](Optics::clustering_order) method
    pub fn reacheability_distance(&self) -> &Vec<f64> { &self.reacheability }

    /// Read-only access to the order in which points were clustered,
    /// provides results of the most recent [`run_clustering()`](Optics::run_clustering) call
    pub fn clustering_order(&self) -> &Vec<usize> { &self.clustering_order }

    /// Copy of clusters created by the most recent [`run_clustering()`](Optics::run_clustering) call.
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
            // --- this will update reacheability distances and insert core points to the list of seeds
            self.update_neighbors(&nb);
            while !self.seeds.is_empty() {
                // --- pop the point with currently lowest reacheability distance value from the priority queue
                let pi = self.seeds[0];
                self.seeds[0] = self.seeds[self.seeds.len() - 1];
                self.seeds.truncate(self.seeds.len() - 1);
                self.seeds.sort_by(|lhs, rhs| self.reacheability[*lhs].partial_cmp(&self.reacheability[*rhs]).unwrap());
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
            if self.reacheability[ci] == self.undefined {
                self.clusters.push(Range { start, end: i });
                start = i;
            }
        }
        self.clusters.push(Range { start, end: self.n });
    }

    fn update_neighbors(&mut self, neighbrs: &Vec<PointAndDistance>) {

        let core_dist = neighbrs[self.min_points - 1].d;
        for n in neighbrs {
            if self.processed[n.idx] {continue}
            let new_reach_dist = n.d.max(core_dist);
            if self.reacheability[n.idx] == self.undefined { // n.idx was not inserted to seed yet
                self.seeds.push(n.idx);
            }
            if self.reacheability[n.idx] > new_reach_dist {
                self.reacheability[n.idx] = new_reach_dist; // set or update its reacheability distance
            }
        }
        self.seeds.sort_by(|lhs, rhs| self.reacheability[*lhs].partial_cmp(&self.reacheability[*rhs]).unwrap());
    }

    /// Allocate internal buffers for a given number of points subjected to clustering
    fn allocate(&mut self) {
        self.reacheability = vec![self.undefined; self.n];   // epsilon * 10 is already the "infinite" distance that defines points that can't be reached
        self.processed = vec![false; self.n];
    }
}


struct NeighborList(Vec<PointAndDistance>);


impl fmt::Display for NeighborList {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.0.iter().fold(Ok(()), |result, pt| {
            result.and_then(|_| write!(f, "{}", pt))
        })
    }
}
