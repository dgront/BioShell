extern crate bioshell_numerical;

use bioshell_numerical::{Vec3, Rototranslation};

#[test]
fn rototranslation_from_column_vectors() {
    let cx: Vec3 = Vec3::new(1.0, 2.0, 3.0);
    let cy: Vec3 = Vec3::new(1.0, 2.0, 3.0);
    let cz: Vec3 = Vec3::new(1.0, 2.0, 3.0);
    let center: Vec3 = Vec3::new(1.0, 2.0, 3.0);

    // let rotot:Rototranslation = Rototranslation::from_column_vectors(cx, cy, cz, center);

    // assert_eq!()
}