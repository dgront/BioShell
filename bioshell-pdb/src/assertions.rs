#[macro_export]
macro_rules! assert_vec3_eq {
    ($va:expr, $vb:expr, $d:expr, $msg:expr) => {
        assert!(($va.x-$vb.x).abs() < $d, "{} : va.x = {}, vb.x = {}", $msg, $va.x, $vb.x);
        assert!(($va.y-$vb.y).abs() < $d, "{} : va.y = {}, vb.y = {}", $msg, $va.y, $vb.y);
        assert!(($va.z-$vb.z).abs() < $d, "{} : va.z = {}, vb.z = {}", $msg, $va.z, $vb.z);
    }
}

#[macro_export]
macro_rules! assert_vec3_ne {
    ($va:expr, $vb:expr, $d:expr, $msg:expr) => {
        assert!(($va.x-$vb.x).abs() > $d, "{} : va.x = {}, vb.x = {}", $msg, $va.x, $vb.x);
        assert!(($va.y-$vb.y).abs() > $d, "{} : va.y = {}, vb.y = {}", $msg, $va.y, $vb.y);
        assert!(($va.z-$vb.z).abs() > $d, "{} : va.z = {}, vb.z = {}", $msg, $va.z, $vb.z);
    }
}

#[macro_export]
macro_rules! assert_delta {
    ($x:expr, $y:expr, $d:expr) => {
        assert!(($x-$y).abs() < $d, "a = {}, b = {}", $x, $y)
    }
}