use std::collections::HashMap;
use std::fmt;

/// Provides one-dimensional histogram with real data counts.
///
/// In this implementation both observed values and their counts are represented as double-precision
/// floating point values. This allows to associate an observation with an arbitrary weight.
///
///  # Examples
///
///  Basic usage:
/// ```
/// use bioshell_statistics::Histogram;
/// let mut h = Histogram::by_bin_width(0.5);
/// h.insert(1.1);
/// h.insert(1.3);
/// assert_eq!(h.which_bin(1.11), 2);
/// assert_eq!(h.get(2), 2.0);
/// ```
#[derive(Debug, Clone)]
pub struct Histogram {
    data: HashMap<i32,f64>,
    bin_width: f64
}

impl Histogram {
    /// Creates a new histogram with a given bin size
    ///
    ///  # Arguments
    /// * `width` - width of each bin of this histogram
    pub fn by_bin_width(width:f64) -> Histogram {
        Histogram{data: HashMap::new(), bin_width:width}
    }

    /// Inserts a value (observation) to this histogram
    pub fn insert(&mut self, x: f64) {
        let bin_id: i32 = (x / self.bin_width).floor() as i32;
        if let Some(x) = self.data.get_mut(&bin_id) {
            *x += 1.0;
        } else {
            self.data.insert(bin_id, 1.0);
        }
    }

    /// Inserts a value (observation) to this histogram
    pub fn insert_weighted(&mut self, x: f64, w: f64) {
        let bin_id: i32 = (x / self.bin_width).floor() as i32;
        if let Some(x) = self.data.get_mut(&bin_id) {
            *x += w;
        } else {
            self.data.insert(bin_id, w);
        }
    }

    /// Says the index of a bin a given value felt into
    pub fn which_bin(&self, val:f64) -> i32 { (val / self.bin_width).floor() as i32 }

    /// Returns the index of the tallest bin
    ///
    /// The maximum value of this histogram may be read by [`get()`] method
    pub fn tallest(&self) -> Option<i32> {
        self.data.iter()
            .max_by(|(_, &value1), (_, &value2)| value1.partial_cmp(&value2).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(&key, _)| key)
    }

    /// Left-hand side range of a given bin
    pub fn bin_min(&self, bin_id: i32) -> f64 { bin_id as f64 * self.bin_width }

    /// Right-hand side range of a given bin (exclusive)
    pub fn bin_max(&self, bin_id: i32) -> f64 { (bin_id + 1) as f64 * self.bin_width }

    /// Returns the count for a given bin
    pub fn get(&self, bin_id: i32) -> f64 { self.data[&bin_id] }

    /// Returns the count for a bin that holds a given value
    pub fn get_by_value(&self, val: f64) -> f64 { self.get(self.which_bin(val)) }

    /// Total number of counts (weights) in this histogram
    pub fn sum(&self) -> f64 { self.data.values().sum() }

    /// Returns the bin width
    pub fn bin_width(&self) -> f64 { self.bin_width }

    /// returns the mode of this histogram
    pub fn mode(&self) -> (f64, f64, f64) {
        let mut max_i: i32 = 0;
        let mut max_v: f64 = 0.0;
        for (i, val) in self.data.iter() {
            if val > &max_v {
                max_i = *i;
                max_v = *val;
            }
        }

        return (self.bin_min(max_i), self.bin_max(max_i), max_v);
    }

    /// Returns observations of this histogram as a vector
    ///
    /// # Example
    /// ```
    /// # use rand::rngs::StdRng;
    /// # use rand::{Rng, SeedableRng};
    /// # use rand_distr::Normal;
    /// # use rand_distr::Distribution;
    /// use bioshell_statistics::Histogram;
    /// let mut rng: StdRng = SeedableRng::from_seed([42; 32]);
    /// let mut normal_distribution = Normal::new(0.0, 1.0).unwrap();
    /// let mut h = Histogram::by_bin_width(0.2);
    /// (0..100).map(|_| normal_distribution.sample(&mut rng)).for_each(|x| h.insert(x));
    /// let values_3sigma = h.to_vector(-3.0, 3.0);
    /// # assert_eq!(h.which_bin(-3.0), -15);
    /// # assert_eq!(h.which_bin(3.0), 15);
    /// assert_eq!(values_3sigma.len(), 31);
    /// ```
    pub fn to_vector(&self, min_val: f64, max_val: f64) -> Vec<f64> {
        let min_key = self.which_bin(min_val);
        let max_key = self.which_bin(max_val);
            (min_key..=max_key)
                .map(|key| self.data.get(&key).cloned().unwrap_or(0.0))
                .collect()
    }

    /// Displays this histogram as an asci-art drawing
    /// ```
    /// # use rand::distributions::Distribution;
    /// # use rand::rngs::StdRng;
    /// # use rand::SeedableRng;
    /// # use rand_distr::Normal;
    /// use bioshell_statistics::Histogram;
    ///
    /// let expected_result =
    /// ".................................................#...................................................
    /// .........................................#################...........................................
    /// .....................................##########################......................................
    /// .................................##################################..................................
    /// ..............................########################################...............................
    /// ..........................###############################################............................
    /// .......................######################################################........................
    /// ...................##############################################################....................
    /// ..............########################################################################...............
    /// .......######################################################################################........
    /// ";
    /// let mut rng: StdRng = SeedableRng::from_seed([42; 32]);
    /// let normal_distribution = Normal::new(0.0, 1.0).unwrap();
    /// let mut h = Histogram::by_bin_width(0.05);
    /// (0..1000000).map(|_| normal_distribution.sample(&mut rng)).for_each(|x| h.insert(x));
    /// let actual = format!("{}", h.draw_horizonaly(-2.5, 2.5, 10));
    /// assert_eq!(actual, expected_result);
    /// ```
    pub fn draw_horizonaly<'a>(&'a self, x_from: f64, x_to: f64, max_height: usize) -> impl fmt::Display + 'a {
        struct Drawing<'a>(&'a Histogram, usize, f64, f64);
        impl<'a> fmt::Display for Drawing<'a> {
            fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
                // Find the maximum bin in the histogram
                let tallest = self.0.tallest().unwrap();
                let max_value = self.0.get(tallest);

                // Calculate the scaling factor for the bars
                let scale_factor = self.1 as f64 / max_value;
                // --- convert histogram to array
                let data = self.0.to_vector(self.2, self.3);
                // Iterate over rows (height) in reverse order to draw vertically
                for h in (0..self.1).rev() {
                    // Iterate over values in the data array
                    for value in &data {
                        // Calculate the height of the bar for the current value
                        let bar_height = (value * scale_factor) as usize;

                        // Print '#' if the current height matches the bar height
                        if bar_height > h {
                            write!(f, "#")?;
                        } else {
                            write!(f, ".")?;
                        }
                    }

                    // Move to the next line for the next row
                    writeln!(f)?;

                }
                Ok(())
            }
        }
        Drawing(self, max_height, x_from, x_to)
    }
}


impl fmt::Display for Histogram {
    /// Creates a `String` representation of a given `Histogram`
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        let mut out: String = String::new();
        for (i, val) in self.data.iter() {
            out = format!("{} {:5}..{:5} [{:4}] {}\n",
                          out, self.bin_min(*i), self.bin_max(*i), i, val);
        }
        write!(f, "{}", out)
    }
}