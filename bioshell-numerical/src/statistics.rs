use std::collections::HashMap;
use std::fmt;

#[derive(Debug, Clone)]
pub struct Histogram {
    data: HashMap<i32,f64>,
    bin_width: f64
}

impl Histogram {
    pub fn by_bin_width(width:f64) -> Histogram {
        Histogram{data: HashMap::new(), bin_width:width}
    }

    pub fn insert(&mut self, x: f64) {
        let bin_id: i32 = (x / self.bin_width).round() as i32;
        if let Some(x) = self.data.get_mut(&bin_id) {
            *x += 1.0;
        } else {
            self.data.insert(bin_id, 1.0);
        }
    }

    pub fn which_bin(&self, val:f64) -> i32 { (val / self.bin_width).round() as i32 }

    pub fn bin_min(&self, bin_id: i32) -> f64 { bin_id as f64 * self.bin_width }

    pub fn bin_max(&self, bin_id: i32) -> f64 { (bin_id + 1) as f64 * self.bin_width }

    pub fn get(&self, bin_id: &i32) -> f64 { self.data[bin_id] }

    pub fn get_by_value(&self, val: f64) -> f64 { self.get(&(self.which_bin(val))) }
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