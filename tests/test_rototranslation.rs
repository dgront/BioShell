extern crate bioshell_numerical;

use bioshell_numerical::{Vec3, Rototranslation};
use bioshell_numerical::RotoTranslation;

#[test]
fn RotoTranslation____from_column_vectors()
{
    let cx: Vec3 = Vec3(1,2,3);
    let cy: Vec3 = Vec3(1,2,3);
    let cz: Vec3 = Vec3(1,2,3);
    let center: Vec3 = Vec3(1,2,3);

    let rotot:Rototranslation = Rototranslation.from_column_vectors(cx, cy, cz, center);

    assert_eq!()
}