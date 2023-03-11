extern crate bioshell_numerical;

use bioshell_numerical::{Vec3, Rototranslation};

fn main()
{
    let _v1: Vec3 = Vec3::new(2.0, 1.0, 3.0);
    let _v2: Vec3 = Vec3::new(1.5, 2.0, 0.3);
    let _v3: Vec3 = Vec3::new(3.7, 1.5, 3.1);
    let _center: Vec3 = Vec3::new(1.0, 1.0, 1.0);
    //let rot1: Rototranslation = Rototranslation::from_column_vectors(&v1, &v2, &v3, &center);
    //rot1.show();
    //let rot2: Rototranslation = Rototranslation::from_row_vectors(&v1, &v2, &v3, &center);
    //rot2.show();
    //let rot3 = Rototranslation::around_axis(&v1, &v2, &v3, 30.0);
    //rot3.show();
    //v1 = v2 - v3;
    //v1.show();
}