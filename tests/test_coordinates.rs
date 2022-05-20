use bioshell_ff::Coordinates;

#[test]
fn coordinates_methods() {
    let mut xyz = Coordinates::new(10);

    assert_eq!(xyz.size(),10);
    xyz[9].x = 9.1;
    assert!(f32::abs(xyz[9].x -9.1) < 1e-5);
}
