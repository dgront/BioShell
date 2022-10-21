use std::cmp::*;
use std::fmt;
use std::ops::{Range};

use crate::clustering::Distance;

/// Provides the OPTICS (ordering points to identify the clustering structure) clustering algorithm.
pub struct Optics<M> {
    /// The radius of a neighborhood.
    pub eps: f64,

    /// The minimum number of points required to form a dense region.
    pub min_samples: usize,

    distance: M,
    n: usize,
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


impl<M: Distance> Optics<M> {
    pub fn new(eps: f64, min_samples: usize, distance: M) -> Self {
        Self { eps, min_samples, distance, n:0, undefined: eps*10.0,
            reacheability: vec![], clustering_order: vec![],
            processed: vec![], seeds: vec![], clusters: vec![] }
    }

    pub fn run_clustering(&mut self, points: &Vec<Vec<f64>>) {
        self.n = points.len();
        self.allocate();

        // --- loop over all points in the dataset
        for i in 0..self.n {
            if self.processed[i] { continue }      // --- skip already processed points
            // --- mark the point as processed, add it to the clustering order
            self.processed[i] = true;
            self.clustering_order.push(i);
            // --- find its neighbors
            let nb = self.find_neighbors(i, points);
            // --- skip the point if it's not a core point
            if nb.len() < self.min_samples {continue}

            // --- ensure the list of seeds is empty (it should be actually)
            self.seeds.clear();
            // --- this will update reacheability distances and insert core points to the list of seeds
            self.update(&nb);
            while !self.seeds.is_empty() {
                // --- pop the point with currently lowest reacheability distance value from the priority queue
                let pi = self.seeds[0];
                self.seeds[0] = self.seeds[self.seeds.len() - 1];
                self.seeds.truncate(self.seeds.len() - 1);
                self.seeds.sort_by(|lhs, rhs| self.reacheability[*lhs].partial_cmp(&self.reacheability[*rhs]).unwrap());
                // --- visit it
                self.clustering_order.push(pi);
                self.processed[pi] = true;
                // --- find its neighbors
                let neighbrs = self.find_neighbors(pi, &points);
                if neighbrs.len() < self.min_samples { continue; }     // --- the pi point is not a core point

                // let nbors:Neighbors = Neighbors{ 0:  neighbrs.clone()};
                // println!("nbors of {}: {}",pi,nbors);

                // --- update the neighbors
                self.update(&neighbrs);
            }
        }

        assert_eq!(self.n, self.clustering_order.len());
        let mut start: usize = 0;
        for i in 1..self.n {
            if self.reacheability[i] == self.undefined {
                self.clusters.push(Range { start, end: i });
                start = i;
            }
        }
        self.clusters.push(Range { start, end: self.n });
    }

    /// Read-only access to the reacheability distance for each point; listed in the order defined
    /// by the [`clustering_order()`](Optics::clustering_order) method
    pub fn reacheability_distance(&self) -> &Vec<f64> { &self.reacheability }

    /// Read-only access to the order in which points were clustered,
    /// provides results of the most recent [`run_clustering()`](Optics::run_clustering) call
    pub fn clustering_order(&self) -> &Vec<usize> { &self.clustering_order }

    /// Read-only access to clusters created by the most recent
    /// [`run_clustering()`](Optics::run_clustering) call
    pub fn clusters(&self) -> &Vec<Range<usize>> { &self.clusters }

    fn update(&mut self, neighbrs: &Vec<PointAndDistance>) {

        let core_dist = neighbrs[self.min_samples - 1].d;
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

    /// Finds epsilon-neighbors of a given point; distances are recorded in the self.dist_storage vector
    fn find_neighbors(&mut self, idx: usize, points: &Vec<Vec<f64>>) -> Vec<PointAndDistance> {
        let mut out: Vec<PointAndDistance> = vec![];
        for i in 0..points.len() {
            let d = self.distance.evaluate(&points[idx], &points[i]);
            if d <= self.eps {
                out.push(PointAndDistance { idx: i, d });
            }
        }
        // --- sort points by distance
        out.sort_by(|k, l| k.d.partial_cmp(&l.d).unwrap());
        return out;
    }

}



#[derive(Clone)]
struct PointAndDistance {
    idx: usize,
    d:f64
}

impl fmt::Display for PointAndDistance {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({}:{:4.3}) ", self.idx, self.d)
    }
}

struct Neighbors(Vec<PointAndDistance>);


impl fmt::Display for Neighbors {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.0.iter().fold(Ok(()), |result, pt| {
            result.and_then(|_| write!(f, "{}", pt))
        })
    }
}
