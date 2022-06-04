use rand::Rng;

use bioshell_ff::Coordinates;

#[test]
fn coordinates_methods() {
    let mut xyz = Coordinates::new(10);

    assert_eq!(xyz.size(), 10);
    xyz[9].x = 9.1;
    assert!(f32::abs(xyz[9].x - 9.1) < 1e-5);
}

#[test]
fn coordinates_pbc() {
    let mut rng = rand::thread_rng();
    const L: f32 = 15.0;

    let mut xyz = Coordinates::new(1);
    xyz.set_box_len(L);
    let l: f32 = xyz.box_len();
    xyz.set(0, rng.gen_range(0.0..l), rng.gen_range(0.0..l), rng.gen_range(0.0..l));

    assert!(xyz.x(0) != xyz.y(0));
    assert!(xyz.x(0) != xyz.z(0));

    for _ in 0..20 {
        assert!(xyz.x(0) >= 0.0);
        assert!(xyz.y(0) >= 0.0);
        assert!(xyz.z(0) >= 0.0);
        assert!(xyz.x(0) <= L);
        assert!(xyz.y(0) <= L);
        assert!(xyz.z(0) <= L);
        xyz.add(0, 0.8, 0.7, 0.9);
    }
}

#[test]
fn closest_distance() {
    let mut rng = rand::thread_rng();
    const L: f32 = 20.0;
    const N: usize = 100;
    const MAX_D: f32 = (L / 2.0) * (L / 2.0) * 3.0;

    let mut xyz = Coordinates::new(N);
    xyz.set_box_len(L);
    let l: f32 = xyz.box_len();

    for i in 0..N {
        xyz.set(i, rng.gen_range(0.0..l), rng.gen_range(0.0..l), rng.gen_range(0.0..l));
        for j in 0..i {
            assert!(xyz.closest_distance_square(i, j) <= MAX_D);
            assert!(xyz.closest_distance_square(i, j) > 0.0);
        }
    }

    xyz.set_box_len(10000000.0);    // --- make all the points well inside the box so closest distance equals to the true distance
    for i in 0..N {
        for j in 0..i {
            assert!((xyz.closest_distance_square(i, j) - xyz.distance_square(i, j)).abs() < 0.0000001);
        }
    }
}
