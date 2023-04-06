/// A set of points subjected to clustering.
///
/// The points themselves must provide a distance metrics that is used during the clustering calculations
pub trait PointsWithDistance {

    /// Returns the distance between two given points.
    ///
    /// # Arguments
    /// * `i` - index of the first point
    /// * `j` - index of the second point
    fn distance(&self, i:usize, j:usize) -> f64;

    /// Returns the total number of data points in this set
    fn count_points(&self) -> usize;
}

impl EuclideanPoints {
    /// Creates a new [`EuclideanPoints`](EuclideanPoints) object from a given vector of points
    /// The input data structure is consumed by this process (i.e. moved)
    ///
    /// # Arguments
    /// * `data` - a 2D vector of f64 values, representing the input data
    ///
    /// # Examples
    /// ```rust
    /// use bioshell_clustering::{EuclideanPoints, PointsWithDistance};
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

impl PointsWithDistance for EuclideanPoints {

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