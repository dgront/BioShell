#[macro_export]
macro_rules! wrap_coordinate_to_box {
    ($val:expr, $L:expr, $coord:expr) => {
        $coord = $val;
        if $coord > $L {
            $coord = $coord - $L
        } else {
            if $coord < 0.0 {
                $coord = $L + $coord
            }
        }
    };
}

/// Calculates the shortest difference `c1 - c2`, taking periodic boundary condition into account
#[macro_export]
macro_rules! closest_image {
    ($c1:expr, $c2:expr, $L: expr,$L2: expr, $delta:expr) => {
        $delta = $c1 - $c2;
        if $delta > 0.0 {
            if $delta > $L2 {
                $delta -= $L
            }
        } else {
            if $delta < -$L2 {
                $delta += $L
            }
        }
    };
}