use bioshell_numerical::statistics::Histogram;

#[test]
fn create_histogram() {
    let test_data = vec!(1.0, 1.1, 1.3, 1.6, 1.7, 2.0);
    let mut h: Histogram = Histogram::by_bin_width(0.5);
    for x in test_data { h.insert(x); }
    assert_eq!(h.which_bin(1.11), 2);
    assert_eq!(h.which_bin(1.49), 2);
    assert_eq!(h.which_bin(1.51), 3);
    assert_eq!(h.sum(), 6.0);
}
