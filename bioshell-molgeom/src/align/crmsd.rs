use nalgebra::{Matrix3, Vector3};
use bioshell_core::{Vec3, Matrix3x3};
use bioshell_pdb::calc::{Rototranslation};

/// Calculate the geometric center of a set of 3D points.
///
/// This is also the center of mass if all points have equal mass.
pub fn centroid(points: &[Vec3]) -> Vec3 {
    let mut c = Vec3::default();
    for p in points { c += p; }
    c /= points.len() as f64;

    return c;
}


/// Calculate the CRMSD between two sets of points.
pub fn crmsd(a: &[Vec3], b: &[Vec3]) -> f64 {
    let (rmsd, _, _, _) = crmsd_core(a, b);
    return rmsd;
}

pub fn crmsd_transform(a: &[Vec3], b: &[Vec3]) -> (f64, Rototranslation) {
    let (rmsd, rotation, centroid_a, centroid_b) = crmsd_core(a, b);

    let mat: [f64; 9] = [
        rotation[(0, 0)], rotation[(0, 1)], rotation[(0, 2)],
        rotation[(1, 0)], rotation[(1, 1)], rotation[(1, 2)],
        rotation[(2, 0)], rotation[(2, 1)], rotation[(2, 2)],
    ];

    let rot = Matrix3x3::from_array(mat);

    let vec_a = Vec3::new(centroid_a.x, centroid_a.y, centroid_a.z);
    let vec_b = Vec3::new(centroid_b.x, centroid_b.y, centroid_b.z);

    let rotated_a = Matrix3x3::mul_vec_s(&rot, &vec_a);

    let mut t = vec_b;
    t -= &rotated_a;

    let rt = Rototranslation::new(rot, t);

    (rmsd, rt)
}

/// Calculate the CRMSD between two sets of points, along with the optimal rotation and centroids.
///
/// Returns a tuple containing the CRMSD, the optimal rotation matrix, and the centroids of both point sets.
fn crmsd_core(a: &[Vec3], b: &[Vec3]) -> (f64, Matrix3<f64>, Vec3, Vec3) {
    assert_eq!(a.len(), b.len(), "Point sets must have the same length");
    assert!(!a.is_empty(), "Point sets must not be empty");

    let n = a.len() as f64;

    let centroid_a = centroid(a);
    let centroid_b = centroid(b);

    let mut h = Matrix3::<f64>::zeros();

    let mut e_a = 0.0;
    let mut e_b = 0.0;

    for (pa, pb) in a.iter().zip(b.iter()) {
        let ax = pa.x - centroid_a.x;
        let ay = pa.y - centroid_a.y;
        let az = pa.z - centroid_a.z;

        let bx = pb.x - centroid_b.x;
        let by = pb.y - centroid_b.y;
        let bz = pb.z - centroid_b.z;

        let va = Vector3::new(ax, ay, az);
        let vb = Vector3::new(bx, by, bz);

        h += va * vb.transpose();

        e_a += ax * ax + ay * ay + az * az;
        e_b += bx * bx + by * by + bz * bz;
    }

    let svd = h.svd(true, true);

    let u = svd.u.expect("SVD failed to compute U");
    let v_t = svd.v_t.expect("SVD failed to compute V^T");

    let v = v_t.transpose();

    let det_sign = if (v * u.transpose()).determinant() < 0.0 {
        -1.0
    } else {
        1.0
    };

    let mut d = Matrix3::<f64>::identity();

    if det_sign < 0.0 {
        d[(2, 2)] = -1.0;
    }

    let rotation = v * d * u.transpose();

    let s = svd.singular_values;
    let optimal_trace = s[0] + s[1] + det_sign * s[2];

    let rmsd_sq = (e_a + e_b - 2.0 * optimal_trace) / n;
    let rmsd = rmsd_sq.max(0.0).sqrt();

    return (rmsd, rotation, centroid_a, centroid_b);
}


