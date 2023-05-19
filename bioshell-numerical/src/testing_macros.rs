/// Equality test for floating point values within a given tolerance
#[macro_export]
macro_rules! assert_eq_float {
    ($lhs:expr, $rhs:expr, $tolerance:expr) => {
        assert!(($lhs-$rhs).abs() < $tolerance, "Floating point comparison between {} and {} failed for tolerance {}", $lhs, $rhs, $tolerance);
    };
}

/// Equality test for Vec3 type within a given tolerance
#[macro_export]
macro_rules! assert_eq_vec3 {
    ($lhs:expr, $rhs:expr, $tolerance:expr) => {
        assert_eq_float!($lhs.x, $rhs.x, $tolerance);
        assert_eq_float!($lhs.y, $rhs.y, $tolerance);
        assert_eq_float!($lhs.z, $rhs.z, $tolerance);
    };
}
