/// Compares two [`Vec3`](crate::calc::Vec3) structs and fails when they differs
///
/// This macro fails an assertion when any of `va` coordinates differs from the corresponding
/// coordinate of `vb` vector by more than `delta`
///
/// # Example
/// ```
/// use bioshell_pdb::assert_vec3_eq;
/// use bioshell_pdb::calc::Vec3;
/// let a = Vec3::new(1.0, 2.0, 3.0);
/// let b = Vec3::new(1.1, 2.0, 3.0);
/// assert_vec3_eq!(a, b, 0.11, "This should not fail");
/// ```
#[macro_export]
macro_rules! assert_vec3_eq {
    ($va:expr, $vb:expr, $delta:expr, $msg:expr) => {
        assert!(($va.x-$vb.x).abs() < $delta, "{} : va.x = {}, vb.x = {}", $msg, $va.x, $vb.x);
        assert!(($va.y-$vb.y).abs() < $delta, "{} : va.y = {}, vb.y = {}", $msg, $va.y, $vb.y);
        assert!(($va.z-$vb.z).abs() < $delta, "{} : va.z = {}, vb.z = {}", $msg, $va.z, $vb.z);
    }
}


/// Compares two floating-point values with tolerance
///
/// This assertion macro fails when the two given values `a` and `b` differ by more than `delta`
#[macro_export]
macro_rules! assert_delta {
    ($a:expr, $b:expr, $delta:expr) => {
        assert!(($a-$b).abs() < $delta, "a = {}, b = {}", $a, $b)
    }
}