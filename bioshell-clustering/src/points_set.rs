use std::fmt;
use std::fmt::Debug;
use std::ops::Index;
use crate::euclidean_distance_squared;
use bioshell_datastructures::kd_tree::{create_kd_tree, find_within, KdTreeNode};

/// Holds information about a neighbor: its index and the distance to it.
///
/// Objects of this class are created by the
/// [`neighbors_of()`](NeighborsOf::neighbors_of) method of the [`NeighborsOf`](NeighborsOf) trait
/// and are used during the OPTICS clustering algorithm.
#[derive(Clone)]
pub struct Neighbor {
    /// index of a neighbor point
    pub idx: usize,
    /// distance to that neighbor
    pub d:f64
}

/// Prints information about a neighbor nicely
impl fmt::Display for Neighbor {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({}:{:4.3}) ", self.idx, self.d)
    }
}

pub trait NeighborsOf {
    /// Returns a vector of neighbors located close enough to a given point
    ///
    /// # Arguments
    /// * `idx` - index of the center point
    /// * `cutoff` - distance cutoff
    /// # Returns
    /// A vector of [`Neighbor`](Neighbor) data structures that provides all points located
    ///no further from `idx` point than the given `cutoff` value
    fn neighbors_of(&self, idx: usize, cutoff: f64) -> Vec<Neighbor>;
}


/// Provides the distance between two points given their indexes.
///
/// In certain clustering scenarios it might be more convenient to provide a struct that returns
/// a distance for two indexes rather than a function that actually computes that. E.g. such an approach
/// simplifies implementation of the OPTICS algorithm.
///
pub trait DistanceByIndex {

    /// Returns the distance between two points given by their indexes.
    ///
    /// # Arguments
    /// * `i` - index of the first point
    /// * `j` - index of the second point
    fn distance(&self, i:usize, j:usize) -> f64;

    /// Returns the total number of data points in this set.
    ///
    /// An index used passed an argument to [`distance()`](distance()) method cannot be greater
    /// that the value returned by this method
    fn count_points(&self) -> usize;
}


pub struct EuclideanPoints<T> {
    dimensionality: usize,
    points: Vec<T>,
    kdtree: Box<KdTreeNode<T>>
}

impl <T>  EuclideanPoints<T> where T: Index<usize, Output = f64> {
    /// Creates a new [`EuclideanPoints`](EuclideanPoints) object from a given vector of points
    /// The input data structure is consumed by this process (i.e. moved)
    ///
    /// # Arguments
    /// * `data` - a 2D vector of f64 values, representing the input data
    ///
    /// # Examples
    /// ```rust
    /// use bioshell_clustering::{EuclideanPoints, DistanceByIndex};
    /// let vecs: Vec<Vec<f64>> = vec![vec![0.0, 0.0], vec![0.5, 1.0], vec![1.5, 0.8]];
    /// let pts = EuclideanPoints::new(vecs, 2);
    /// assert!((1.118-pts.distance(0, 1)).abs() < 0.001);
    ///
    /// let arrays: Vec<[f64; 2]> = vec![[0.0, 0.0], [0.5, 1.0], [1.5, 0.8]];
    /// let pts = EuclideanPoints::new(arrays, 2);
    /// assert!((1.118-d.distance(0, 1)).abs() < 0.001);
    /// ```
    pub fn new(data: Vec<T>, dimensionality: usize) -> EuclideanPoints<T> where T: Clone+Debug {
        let mut data_copy = data.clone();
        let root = create_kd_tree(&mut data_copy, dimensionality);
        EuclideanPoints {dimensionality, points: data, kdtree: root.unwrap() }
    }
}

impl <T> DistanceByIndex for EuclideanPoints<T> where T: Index<usize, Output = f64> {

    fn distance(&self, i: usize, j: usize) -> f64 {
        euclidean_distance_squared(&self.points[i], &self.points[j], self.dimensionality)
    }

    fn count_points(&self) -> usize { self.points.len() }
}

impl<T> NeighborsOf for EuclideanPoints<T> where T: Index<usize, Output = f64> {

    fn neighbors_of(&self, idx:usize, cutoff: f64) -> Vec<Neighbor> {

        let res = find_within(&self.kdtree, &self.points[idx], self.dimensionality,
                    cutoff*cutoff, euclidean_distance_squared);
        let mut out: Vec<Neighbor> = vec![];
        for i in 0..self.count_points() {
            let d = self.distance(idx, i);
            if d <= cutoff {
                out.push(Neighbor { idx: i, d });
            }
        }
        // --- sort points by distance
        out.sort_by(|k, l| k.d.partial_cmp(&l.d).unwrap());
        return out;
    }
}

/*
impl EuclideanPoints {
    /// Creates a new [`EuclideanPoints`](EuclideanPoints) object from a given vector of points
    /// The input data structure is consumed by this process (i.e. moved)
    ///
    /// # Arguments
    /// * `data` - a 2D vector of f64 values, representing the input data
    ///
    /// # Examples
    /// ```rust
    /// use bioshell_clustering::{EuclideanPoints, DistanceByIndex};
    /// let points: Vec<Vec<f64>> = vec![vec![0.0, 0.0], vec![0.5, 1.0], vec![1.5, 0.8]];
    /// let d = EuclideanPoints::new(points);
    /// assert!((1.118-d.distance(0, 1)).abs() < 0.001);
    /// ```
    pub fn new(data: Vec<Vec<f64>>) -> EuclideanPoints {
        EuclideanPoints {datapoints: data}
    }
}


/// A container for N-dimensional points of `Vec<f64>` type and Euclidean distance
pub struct EuclideanPoints {
    datapoints: Vec<Vec<f64>>
}

impl DistanceByIndex for EuclideanPoints {

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

    fn count_points(&self) -> usize { self.datapoints.len() }
}

 */