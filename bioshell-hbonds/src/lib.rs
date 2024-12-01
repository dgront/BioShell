//! Discovers and describes hydrogen bonds in protein molecules.

mod backbone_hb_map;
pub use backbone_hb_map::*;

/// Define the maximum allowed distance between a hydrogen atom and its acceptor to still record a hydrogen bond
const MAX_AH_DISTANCE: f64 = 3.0;

/// Define the maximum allowed distance between the two heavy atoms: donor and acceptor to still record a hydrogen bond
const MAX_AD_DISTANCE: f64 = 4.0;

/// The minimum value of O..H-N angle to exist in a hydrogen bond (in degrees). If the angle is smaller, a H-bond is not detected
const MIN_AHD_ANGLE: f64 = 100.0;

/// The minimum value of C=O..H angle to exist in a hydrogen bond (in degrees). If the angle is smaller, a H-bond is not detected
const MIN_PAH_ANGLE: f64 = 100.0;

/// Length of the amide N-H bond, used in reconstruction of amide H position
const N_H_BOND_LENGTH: f64 = 1.0;