use nalgebra::{Matrix3, Vector3};
use bioshell_pdb::calc::Vec3;

pub fn centroid(points: &[Vec3]) -> Vec3 {
    let mut c = Vec3::default();
    for p in points { c += p; }
    c /= points.len() as f64;

    return c;
}


pub fn crmsd(a: &[Vec3], b: &[Vec3]) -> f64 {

    assert_eq!(a.len(), b.len(), "Point sets must have the same length");
    assert!(!a.is_empty(), "Point sets must not be empty");

    let n = a.len() as f64;

    let centroid_a = centroid(a);
    let centroid_b = centroid(b);

    let mut h = Matrix3::<f64>::zeros();

    for (pa, pb) in a.iter().zip(b.iter()) {
        let mut va = *pa;
        let mut vb = *pb;

        va -= &centroid_a;
        vb -= &centroid_b;

        let va = Vector3::new(va.x, va.y, va.z);
        let vb = Vector3::new(vb.x, vb.y, vb.z);

        h += va * vb.transpose();
    }

    let svd = h.svd(true, true);

    let u = svd.u.expect("SVD failed to compute U");
    let v_t = svd.v_t.expect("SVD failed to compute V^T");

    let mut d = Matrix3::<f64>::identity();

    if (v_t.transpose() * u.transpose()).determinant() < 0.0 {
        d[(2, 2)] = -1.0;
    }

    let rotation = v_t.transpose() * d * u.transpose();

    let mut sum_sq = 0.0;

    for (pa, pb) in a.iter().zip(b.iter()) {
        let mut va = *pa;
        let mut vb = *pb;

        va -= &centroid_a;
        vb -= &centroid_b;

        let va = Vector3::new(va.x, va.y, va.z);
        let vb = Vector3::new(vb.x, vb.y, vb.z);

        let diff = rotation * va - vb;
        sum_sq += diff.norm_squared();
    }

    return (sum_sq / n).sqrt();
}