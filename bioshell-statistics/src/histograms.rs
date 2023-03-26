use std::collections::HashMap;
use std::fmt;

/// Provides one-dimensional histogram with real data counts
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
    data: HashMap<i32, f64>,
    bin_width: f64,
}

impl Histogram {
    /// Creates a new histogram with a given bin size
    ///
    ///  # Arguments
    /// * `width` - width of each bin of this histogram
    pub fn by_bin_width(width: f64) -> Histogram {
        Histogram {
            data: HashMap::new(),
            bin_width: width,
        }
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

    /// Says the index of a bin a given value felt into
    pub fn which_bin(&self, val: f64) -> i32 {
        (val / self.bin_width).floor() as i32
    }

    /// Left-hand side range of a given bin
    pub fn bin_min(&self, bin_id: i32) -> f64 {
        bin_id as f64 * self.bin_width
    }

    /// Right-hand side range of a given bin (exclusive)
    pub fn bin_max(&self, bin_id: i32) -> f64 {
        (bin_id + 1) as f64 * self.bin_width
    }

    /// Returns the count for a given bin
    pub fn get(&self, bin_id: i32) -> f64 {
        self.data[&bin_id]
    }

    /// Returns the count for a bin that holds a given value
    pub fn get_by_value(&self, val: f64) -> f64 {
        self.get(self.which_bin(val))
    }

    /// Total number of counts in this histogram
    pub fn sum(&self) -> f64 {
        self.data.values().sum()
    }

    /// Returns the bin width
    pub fn bin_width(&self) -> f64 {
        self.bin_width
    }

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
}

impl fmt::Display for Histogram {
    /// Creates a `String` representation of a given `Histogram`
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        let mut out: String = String::new();
        for (i, val) in self.data.iter() {
            out = format!("{} {:5}..{:5} [{:4}] {}\n", out, self.bin_min(*i), self.bin_max(*i), i, val);
        }
        write!(f, "{}", out)
    }
}
