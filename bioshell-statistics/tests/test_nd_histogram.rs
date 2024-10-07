#[cfg(test)]
mod tests {
    use bioshell_statistics::HistogramND;
    use super::*;

    #[test]
    fn test_insert() {
        let mut hist = HistogramND::<2>::by_bin_widths([1.0, 0.5]);

        // Insert values and check the bin counts
        hist.insert([0.5, 0.25]);
        assert_eq!(hist.get(&vec![0, 0]), 1.0);

        hist.insert([1.5, 0.75]);
        assert_eq!(hist.get(&vec![1, 1]), 1.0);

        hist.insert([1.5, 0.75]);
        assert_eq!(hist.get(&vec![1, 1]), 2.0);
    }

    #[test]
    fn test_insert_weighted() {
        let mut hist = HistogramND::<3>::by_bin_widths([2.0, 0.5, 1.0]);

        // Insert weighted values and check the bin counts
        hist.insert_weighted([3.5, 0.75, 2.0], 1.5);
        assert_eq!(hist.get(&vec![1, 1, 2]), 1.5);

        hist.insert_weighted([3.5, 0.75, 2.0], 0.5);
        assert_eq!(hist.get(&vec![1, 1, 2]), 2.0);
    }

    #[test]
    fn test_which_bin() {
        let hist = HistogramND::<2>::by_bin_widths([1.0, 0.5]);

        // Check the bin indices for values
        assert_eq!(hist.which_bin(&[0.5, 0.25]), vec![0, 0]);
        assert_eq!(hist.which_bin(&[1.5, 0.75]), vec![1, 1]);
        assert_eq!(hist.which_bin(&[2.9, -0.5]), vec![2, -1]);
    }

    #[test]
    fn test_mode() {
        let mut hist = HistogramND::<2>::by_bin_widths([1.0, 0.5]);

        // Insert values into different bins
        hist.insert([0.5, 0.25]);  // Bin [0, 0]
        hist.insert([1.5, 0.75]);  // Bin [1, 1]
        hist.insert([1.5, 0.75]);  // Bin [1, 1]
        hist.insert([2.5, -0.25]); // Bin [2, -1]

        // Get the mode (tallest bin)
        let (min, max, count) = hist.mode();
        assert_eq!(min, vec![1.0, 0.5]);
        assert_eq!(max, vec![2.0, 1.0]);
        assert_eq!(count, 2.0);
    }

    #[test]
    fn test_sum() {
        let mut hist = HistogramND::<2>::by_bin_widths([1.0, 0.5]);

        // Insert values into different bins
        hist.insert([0.5, 0.25]);
        hist.insert([1.5, 0.75]);
        hist.insert([2.5, -0.25]);

        // Check the sum of all bins
        assert_eq!(hist.sum(), 3.0);
    }


    #[test]
    fn test_max() {
        let mut hist = HistogramND::<2>::by_bin_widths([1.0, 0.5]);

        // Insert values into different bins
        hist.insert([0.75, 0.25]);  // Bin [0, 0]
        hist.insert([0.5, 0.25]);  // Bin [0, 0]
        hist.insert([1.5, 0.75]);  // Bin [1, 1]
        hist.insert([2.5, -0.25]); // Bin [2, -1]

        // Check min and max bin indices
        assert_eq!(hist.max(), Some(vec![0, 0]));

        let mut hist = HistogramND::<2>::by_bin_widths([1.0, 0.5]);

        // Insert values into different bins
        hist.insert([0.5, 0.25]);  // Bin [0, 0]
        hist.insert([1.5, 0.75]);  // Bin [1, 1]
        hist.insert([1.5, 0.75]);  // Bin [1, 1]
        hist.insert([2.5, -0.25]); // Bin [2, -1]

        // Check that the tallest bin is the one with the most inserts
        assert_eq!(hist.max(), Some(vec![1, 1]));
    }
}
