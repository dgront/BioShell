use bioshell_core::calc::structure::*;
use bioshell_core::structure::Coordinates;

#[test]
fn test_cm_rg_rend() {
    // ---------- prepare coordinates like that:
    //   _   _
    //  | |_| |
    //  ______|
    // CM should be in (1.5, 0, 0)
    let xy = vec![
        (0f64, 0f64),
        (0.0, 1.0),
        (1.0, 1.0),
        (1.0, 0.0),
        (2.0, 0.0),
        (2.0, 1.0),
        (3.0, 1.0),
        (3.0, 0.0),
        (3.0, -1.0),
        (2.0, -1.0),
        (1.0, -1.0),
        (0.0, -1.0),
    ];
    let mut cords: Coordinates = Coordinates::new(xy.len());

    for i in 0..xy.len() {
        cords[i].x = xy[i].0;
        cords[i].y = xy[i].1;
    }
    // ---------- Calculate CM
    let (x, y, z) = cm(&cords, 0);

    assert_eq!(x, 1.5);
    assert_eq!(y, 0.0);
    assert_eq!(z, 0.0);

    assert_eq!(r_end_squared(&cords, 0), 1.0);
}
