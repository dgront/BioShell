use nalgebra::{DMatrix, SVD, Vector3};
use crate::calc::Vec3;
use crate::HasCartesians;

/// A 3D vector defined of a molecular fragment, e.g. a helical axis.
///
///
#[derive(Debug, Clone)]
pub struct SubstructureAxis {
    versor: Vec3,   // Unit vector along the axis (direction)
    begin: Vec3,    // N-terminal end of the vector
    end: Vec3,      // C-terminal end of the vector
}

impl SubstructureAxis {
    /// Computes the helical axis from a list of ordered positions in 3D space.
    pub fn from_3d_points<T: HasCartesians>(coords: &[T]) -> Self {
        let n = coords.len();
        assert!(n >= 2, "At least two coordinates are required");

        // Internal conversion
        fn to_na(v: &Vec3) -> Vector3<f64> { Vector3::new(v[0], v[1], v[2]) }

        fn from_na(v: &Vector3<f64>) -> Vec3 { Vec3::from_array(&[v[0], v[1], v[2]]) }

        // Step 1: Compute centroid
        let mut centroid = Vec3::from_float(0.0);
        for v in coords {
            centroid += v.position();
        }
        centroid /= n as f64;
        let centroid_na = to_na(&centroid);

        // Step 2: Centered matrix for PCA
        let mut centered_data: Vec<f64> = Vec::with_capacity(3 * n);
        for v in coords {
            let c = Vec3::sub_s(v.position(), &centroid);
            centered_data.push(c[0]);
            centered_data.push(c[1]);
            centered_data.push(c[2]);
        }
        let centered = DMatrix::from_row_slice(n, 3, &centered_data);

        // Step 3: Covariance matrix & SVD
        let cov = &centered.transpose() * &centered / (n as f64);
        let svd = SVD::new(cov, true, false);
        let u = svd.u.expect("SVD failed");
        let axis_dvec = u.column(0); // Principal component
        let mut axis = Vector3::new(axis_dvec[0], axis_dvec[1], axis_dvec[2]);

        // Step 4: Orient from N- to C-terminal
        let end_vec = Vec3::sub_s(&coords[n - 1].position(), &coords[0].position());
        let axis_vec3 = from_na(&axis);
        if Vec3::dot(&axis_vec3, &end_vec) < 0.0 {
            axis *= -1.0;
        }

        let versor_na = axis.normalize();
        let versor = from_na(&versor_na);

        // Step 5: Project all points onto the axis (in scalar form)
        let projections: Vec<f64> = coords.iter()
            .map(|v| (to_na(v.position()) - centroid_na).dot(&versor_na))
            .collect();

        let min_proj = projections.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_proj = projections.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

        // Step 6: Compute true begin and end points (center of circular bases)
        let begin_na = centroid_na + versor_na * min_proj;
        let end_na = centroid_na + versor_na * max_proj;

        Self {
            versor,
            begin: from_na(&begin_na),
            end: from_na(&end_na),
        }
    }

    /// Returns the unit vector along the axis of the helix.
    pub fn versor(&self) -> Vec3 {
        self.versor
    }

    /// Returns the center of the base at the N-terminal end.
    pub fn begin(&self) -> Vec3 {
        self.begin
    }

    /// Returns the center of the base at the C-terminal end.
    pub fn end(&self) -> Vec3 {
        self.end
    }

    pub fn length(&self) -> f64 {
        self.begin.distance_to(&self.end)
    }

    /// Returns the midpoint of the helix axis (center of the cylinder).
    pub fn midpoint(&self) -> Vec3 {
        Vec3::from_array(&[
            0.5 * (self.begin[0] + self.end[0]),
            0.5 * (self.begin[1] + self.end[1]),
            0.5 * (self.begin[2] + self.end[2]),
        ])
    }
}

