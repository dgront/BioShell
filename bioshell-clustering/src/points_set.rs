use std::fmt;
use std::fmt::Debug;
use std::ops::Index;
use bioshell_datastructures::BinaryTreeNode;

use bioshell_datastructures::kd_tree::{create_kd_tree, find_within, KdTreeData};

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

/// Lists neighbors of a given point within a set.
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

#[derive(Clone, Debug)]
struct IndexedPoint<T> {
    index: usize,
    point: T
}

/// Indexing operator delegates to the operator of the inner object of a generic type T
impl<T> Index<usize> for IndexedPoint<T> where T: Index<usize, Output = f64> {
    type Output = f64;

    fn index(&self, index: usize) -> &f64 { &self.point[index] }
}


/// Creates a set of points that will be subjected to clustering analysis.
///
/// [`CartesianPoints`](CartesianPoints) struct holds objects of generic type `T`, which must
/// be indexable (i.e. must implement `Index<usize, Output = f64>` trait. Examples of the type `T`
/// are `[f64;3]`, `Vec<f64>` and `bioshell-numerical::Vec3`. The second generic argument `D`
/// defines the distance function between points `T`, which must be restricted to `Fn(&T, &T, usize) -> f64`
///
/// The [`CartesianPoints`](CartesianPoints) struct implements [`NeighborsOf`](NeighborsOf)
/// and [`DistanceByIndex`](DistanceByIndex) traits which is required for the [`Optics`](Optics) clustering method.
/// [`bioshell-datastructures::kd_tree`](bioshell_datastructures::kd_tree) module is used to
/// efficiently compile a list of  neighbors for each query point
pub struct CartesianPoints<T, D> {
    distance_function: D,
    dimensionality: usize,
    points: Vec<IndexedPoint<T>>,
    kdtree: Box<BinaryTreeNode<KdTreeData<IndexedPoint<T>>>>,
}

impl <T, D>  CartesianPoints<T, D> where T: Index<usize, Output = f64> {
    /// Creates a new [`CartesianPoints`](CartesianPoints) object from a given vector of points
    /// The input data structure is consumed by this process (i.e. moved)
    ///
    /// # Arguments
    /// * `dist_func` - a callable object of type `D` that can be used to compute the distance between two
    ///     objects of type `T`
    /// * `data` - a vector of  objects of type `T`, representing the input data
    /// * `dimensionality` - number of elements (or dimensions) of each object of type `T`
    ///
    /// # Generic types
    /// * `T` - the type of data subjected to clustering; it must be indexable, i.e. it must be possible to
    ///     reach i-th element of `T` with an indexing operator `T[i]` which must return `f64`
    /// * `D` - the type of a distance function, which must be compatible with `Fn(&T, &T, usize) -> f64`;
    ///     it accepts three arguments: two points of type `T` and the number of dimension for the `T` objects
    ///
    ///
    /// # Examples
    /// ```rust
    /// use bioshell_clustering::{CartesianPoints, DistanceByIndex, euclidean_distance_squared};
    /// let vecs: Vec<Vec<f64>> = vec![vec![0.0, 0.0], vec![0.5, 1.0], vec![1.5, 0.8]];
    /// let pts = CartesianPoints::new(euclidean_distance_squared, vecs, 2);
    /// assert!((1.25-pts.distance(0, 1)).abs() < 0.001);
    ///
    /// let arrays: Vec<[f64; 2]> = vec![[0.0, 0.0], [0.5, 1.0], [1.5, 0.8]];
    /// let pts = CartesianPoints::new(euclidean_distance_squared, arrays, 2);
    /// assert!((1.25-pts.distance(0, 1)).abs() < 0.001);
    /// ```
    pub fn new(dist_func: D, data: Vec<T>, dimensionality: usize) -> CartesianPoints<T, D> where T: Clone+Debug {
        let mut ipoints: Vec<IndexedPoint<T>> = vec![];
        for (i, d) in data.iter().enumerate() { ipoints.push(IndexedPoint{ index: i, point: d.clone() }); }
        let mut ipoints_copy = ipoints.clone();
        let root = create_kd_tree(&mut ipoints_copy, dimensionality);
        CartesianPoints {distance_function: dist_func, dimensionality, points: ipoints, kdtree: root.unwrap() }
    }

    /// Returns the dimensionality of points
    pub fn dimensionality(&self) -> usize { self.dimensionality }
}

impl <T, D> DistanceByIndex for CartesianPoints<T, D>
    where T: Index<usize, Output = f64>, D: Fn(&T, &T, usize) -> f64  {

    fn distance(&self, i: usize, j: usize) -> f64 {
        (self.distance_function)(&self.points[i].point, &self.points[j].point, self.dimensionality)
    }

    fn count_points(&self) -> usize { self.points.len() }
}

impl<T, D> NeighborsOf for CartesianPoints<T, D>
    where T: Index<usize, Output = f64>, D: Fn(&T, &T, usize) -> f64  {

    fn neighbors_of(&self, idx:usize, cutoff: f64) -> Vec<Neighbor> {


        let mut out: Vec<Neighbor> = vec![];
        let df = |a: &IndexedPoint<T>, b: &IndexedPoint<T>, dim: usize| {
            (self.distance_function)(&a.point, &b.point, dim)
        };
        let nb = find_within(&self.kdtree, &self.points[idx], self.dimensionality,
                             cutoff, &df);

        for pi in nb {
            let d = self.distance(idx, pi.index);
            if d <= cutoff {
                out.push(Neighbor { idx: pi.index, d });
            }
        }
        // --- sort points by distance
        out.sort_by(|k, l| k.d.partial_cmp(&l.d).unwrap());
        return out;
    }
}
