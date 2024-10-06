use std::collections::HashMap;
use std::ops::Range;

/// Provides N-dimensional histogram of N-dimensional data.
///
/// # Example
/// ```rust
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// # use rand::distributions::Distribution;
/// # use rand::rngs::StdRng;
/// # use rand::SeedableRng;
/// # use rand_distr::Normal;
/// use bioshell_statistics::HistogramND;
///
/// let mut rng: StdRng = SeedableRng::from_seed([42; 32]);
/// let normal = Normal::new(0.0, 1.0)?;
///
/// let mut histogram = HistogramND::<2>::by_bin_widths([0.1, 0.1]);
///
/// for _ in 0..100_000 {
///     let x = normal.sample(&mut rng);
///     let y = normal.sample(&mut rng);
///     histogram.insert([x, y]);
/// }
///
/// // Verify histogram mode (tallest bin)
/// let (min, max, count) = histogram.mode();
/// println!("# Histogram mode (x, y): min = {:?}, max = {:?}, count = {:.4}", min, max, count);
/// # Ok(())
/// # }
/// ```
pub struct HistogramND<const N: usize> {
    data: HashMap<Vec<i32>, f64>,
    bin_widths: [f64; N],
}

impl<const N: usize> HistogramND<N> {
    /// Creates a new N-dimensional histogram with given bin sizes for each dimension
    ///
    ///  # Arguments
    /// * `widths` - array of bin widths for each dimension
    pub fn by_bin_widths(widths: [f64; N]) -> HistogramND<N> {
        HistogramND {
            data: HashMap::new(),
            bin_widths: widths,
        }
    }

    /// Inserts a new observation to this histogram
    ///
    /// # Examples
    /// ```
    /// use bioshell_statistics::HistogramND;
    /// let mut h2 = HistogramND::<2>::by_bin_widths([0.1, 0.1]);
    /// h2.insert([0.1, 0.1]);
    /// let mut h5 = HistogramND::<5>::by_bin_widths([0.1, 0.1, 0.1, 0.1, 0.1]);
    /// h5.insert([0.1, 0.3, 0.5, 0.7, 0.9]);
    /// ```
    pub fn insert(&mut self, x: [f64; N]) {
        let bin_id = self.which_bin(&x);
        if let Some(count) = self.data.get_mut(&bin_id) {
            *count += 1.0;
        } else {
            self.data.insert(bin_id, 1.0);
        }
    }

    /// Inserts a value (observation) with a weight
    pub fn insert_weighted(&mut self, x: [f64; N], w: f64) {
        let bin_id = self.which_bin(&x);
        if let Some(count) = self.data.get_mut(&bin_id) {
            *count += w;
        } else {
            self.data.insert(bin_id, w);
        }
    }

    /// Says the index of a bin a given value fell into
    pub fn which_bin(&self, val: &[f64; N]) -> Vec<i32> {
        val.iter()
            .enumerate()
            .map(|(i, &x)| (x / self.bin_widths[i]).floor() as i32)
            .collect()
    }

    /// Returns the index of a bin that holds the largest value observed by this histogram
    pub fn max(&self) -> Option<Vec<i32>> {
        self.data
            .iter()
            .max_by(|(_, &v1), (_, &v2)| v1.partial_cmp(&v2).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(key, _)| key.clone())
    }

    /// Left-hand side range of a given bin in a specified dimension
    pub fn bin_min(&self, bin_id: &Vec<i32>, dim: usize) -> f64 {
        bin_id[dim] as f64 * self.bin_widths[dim]
    }

    /// Right-hand side range of a given bin (exclusive) in a specified dimension
    pub fn bin_max(&self, bin_id: &Vec<i32>, dim: usize) -> f64 {
        (bin_id[dim] + 1) as f64 * self.bin_widths[dim]
    }

    /// Returns the count for a given bin
    pub fn get(&self, bin_id: &Vec<i32>) -> f64 {
        *self.data.get(bin_id).unwrap_or(&0.0)
    }

    /// Returns the count for a bin that holds a given value
    pub fn get_by_value(&self, val: [f64; N]) -> f64 {
        self.get(&self.which_bin(&val))
    }

    /// Total number of counts (weights) in this histogram
    pub fn sum(&self) -> f64 {
        self.data.values().sum()
    }

    /// Returns the bin widths for each dimension
    pub fn bin_widths(&self) -> &[f64; N] {
        &self.bin_widths
    }

    /// Returns the mode of this histogram
    pub fn mode(&self) -> (Vec<f64>, Vec<f64>, f64) {
        let mut max_bin: Vec<i32> = vec![0; N];
        let mut max_value: f64 = 0.0;

        for (bin, &value) in &self.data {
            if value > max_value {
                max_bin = bin.clone();
                max_value = value;
            }
        }

        let bin_min = (0..N)
            .map(|i| self.bin_min(&max_bin, i))
            .collect::<Vec<f64>>();
        let bin_max = (0..N)
            .map(|i| self.bin_max(&max_bin, i))
            .collect::<Vec<f64>>();

        (bin_min, bin_max, max_value)
    }
}

use std::fmt;

impl fmt::Display for HistogramND<2> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Loop over each bin in the histogram
        for (bin, &count) in &self.data {
            if count > 0.0 {
                // Calculate the min values for the x and y dimensions
                let x_min = self.bin_min(bin, 0);
                let y_min = self.bin_min(bin, 1);
                writeln!(f, "{:.4}\t{:.4}\t{:.4}", x_min, y_min, count)?;
            }
        }
        Ok(())
    }
}

