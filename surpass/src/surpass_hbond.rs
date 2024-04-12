use bioshell_pdb::calc::{Rototranslation, Vec3};
use crate::{MoveProposal, SurpassAlphaSystem, SurpassEnergy};

/// Shortest CA-CA distance the residues might form a hydrogen bond
const MIN_HBOND_DISTANCE_SQ: f64 = 4.0 * 4.0;
/// Longest CA-CA distance the residues might form a hydrogen bond
const MAX_HBOND_DISTANCE_SQ: f64 = 6.2 * 6.2;

pub struct HBond3CA {
      h_along: HelicalHBond,
      h_against: HelicalHBond,
      e_parallel: StrandHBond,
      e_anti: StrandHBond,
}

impl HBond3CA {
    pub fn new() -> HBond3CA {
         let h_along = HelicalHBond::new(-0.2, -1.3,
                                         -0.778, 0.236,
                                         4.492, 0.095,
                                         3.138, 0.069);
        let h_against = HelicalHBond::new(-0.2, -1.3,
                                          -0.778, 0.236,
                                          4.492, 0.095,
                                          3.138, 0.069);
        let e_parallel = StrandHBond::new(-0.2, -1.3,
                                          -0.0747412, 0.300409,
                                          4.79221, 0.16418,
                                          3.2715, 0.175267,
                                          0.0783662, 0.240884,
                                          4.79591, 0.145513,
                                          6.42867, 0.160225,                                                    );
        let e_anti = StrandHBond::new(-0.2, -1.3,
                                      -0.715211, 0.304241,
                                      4.81181, 0.127693,
                                      3.37516, 0.243099,
                                      -0.512375, 0.299921,
                                      4.74849, 0.173994,
                                      6.52424, 0.216474);

        return HBond3CA{ h_along, h_against, e_parallel, e_anti };
    }
    /// Calculate coarse-grained hydrogen bond energy between two residues of SURPASS-alpha system
    ///
    pub fn evaluate_hbond_energy(&self, system: &SurpassAlphaSystem, i_res: usize, j_res: usize) -> f64 {

        // ---------- First and last CA of a chain can't form an H-bond
        let chain_i_range = system.chain_atoms(system.chain(i_res) as usize);
        if i_res == chain_i_range.start || i_res == chain_i_range.end -1 { return 0.0; }
        let chain_j_range = system.chain_atoms(system.chain(j_res) as usize);
        if j_res == chain_j_range.start || j_res == chain_j_range.end -1 { return 0.0; }

        // ---------- get Cartesian coordinates of the 6 alpha carbons correctly wrapped in PBC
        let vi = system.ca_to_vec3(i_res);
        let vj = system.ca_to_vec3(j_res);

        let vip = system.ca_to_nearest_vec3(i_res - 1, i_res);
        let vin = system.ca_to_nearest_vec3(i_res + 1, i_res);
        let vjp = system.ca_to_nearest_vec3(j_res - 1, j_res);
        let vjn = system.ca_to_nearest_vec3(j_res + 1, j_res);

        let out = self.evaluate_hbond_energy_6ca(i_res, &vip, &vi, &vin, j_res, &vjp, &vj, &vjn);

        // println!("{} {}  {}",i_res, j_res, out);
        return out;
    }

    /// Calculate coarse-grained hydrogen bond energy between 6 amino acid residues
    ///
    pub fn evaluate_hbond_energy_6ca(&self, i_res: usize, vip: &Vec3, vi: &Vec3, vin: &Vec3,
                                     j_res: usize, vjp: &Vec3, vj: &Vec3, vjn: &Vec3) -> f64 {

        // ---------- check the sense of a possible  H-bond
        let sense = sense_h_bond(&vip, &vi, &vin, &vjp, &vj, &vjn);
        match sense {
            HBondSense::None => { return 0.0;}
            HBondSense::Parallel => {
                if j_res > i_res && j_res - i_res == 4 {
                    return -(self.h_against.evaluate(&vjp, &vj, &vjn, &vi)*
                        self.h_along.evaluate(&vip, &vi, &vin, &vj)).powf(1.0 / 6.0);
                }
                else if j_res > i_res && j_res - i_res == 4 {
                    return -(self.h_along.evaluate(&vjp, &vj, &vjn, &vi)*
                        self.h_against.evaluate(&vip, &vi, &vin, &vj)).powf(1.0 / 6.0);
                }

                // println!("parallel E");
                return - (self.e_parallel.evaluate(&vip, &vi, &vin, &vj)
                    * self.e_parallel.evaluate(&vjp, &vj, &vjn, &vi)).powf(1.0 / 6.0);
            }
            HBondSense::Antiparallel => {
                // println!("antiparallel E");

                return -(self.e_anti.evaluate(&vip, &vi, &vin, &vj)
                    * self.e_anti.evaluate(&vjn, &vj, &vjp, &vi)).powf(1.0 / 6.0);
            }
        }
    }
}

impl SurpassEnergy for HBond3CA {
    fn evaluate(&self, conf: &SurpassAlphaSystem) -> f64 {
        let mut total_en = 0.0;

        // ---------- Energy of the chain 0 with itself
        let i_range = conf.chain_atoms(0);
        for i_pos in i_range.start+4..i_range.end-1 {
            for j_pos in i_range.start+1..i_pos-2 {
                total_en += self.evaluate_hbond_energy(conf, i_pos, j_pos);
            }
        }

        // ---------- Loop over chain "i" starting from chain 1
        for i_chain in 1..conf.count_chains() {
            let i_range = conf.chain_atoms(i_chain);
            // ---------- Energy of i_chain with itself
            for i_pos in i_range.start+4..i_range.end-1 {
                for j_pos in i_range.start+1..i_pos-2 {
                    total_en += self.evaluate_hbond_energy(conf, i_pos, j_pos);
                }
            }
            // ---------- Energy of i_chain with j_chain, i_chain > j_chain
            for j_chain in 0..i_chain {
                let j_range = conf.chain_atoms(j_chain);
                for i_pos in i_range.start+1..i_range.end-1 {
                    for j_pos in j_range.start+1..j_range.end-1 {
                        total_en += self.evaluate_hbond_energy(conf, i_pos, j_pos);
                    }
                }
            }
        }

        return total_en;
    }

    fn evaluate_delta(&self, conf: &SurpassAlphaSystem, move_prop: &MoveProposal) -> f64 {

        let mut moved_system = conf.clone();
        move_prop.apply(&mut moved_system);
        return self.evaluate(&moved_system) - self.evaluate(&conf);
    }
}

struct HBondFunction {
    x_from: f64,
    x_to: f64,
    r_from: f64,
    r_to: f64,
    a_from: f64,
    a_to: f64,
    smooth: f64,
}

impl HBondFunction {
    pub fn evaluate(&self, x: f64, r: f64, a: f64) -> f64 {
        // ---------- calculate the energy value
        let en_x = ((x - self.x_from) / self.smooth).tanh() - ((x - self.x_to) / self.smooth).tanh();
        let en_r = ((r - self.r_from) / self.smooth).tanh() - ((r - self.r_to) / self.smooth).tanh();
        let en_a = ((a - self.a_from) / self.smooth).tanh() - ((a - self.a_to) / self.smooth).tanh();
        return en_x * en_r * en_a / 8.0;
    }
}
struct HelicalHBond {
    y_center: f64,
    z_center: f64,
    function: HBondFunction
}

impl HelicalHBond {
    pub fn new(y_center: f64, z_center: f64, x_avg: f64, x_std: f64, r_avg: f64, r_std: f64, a_avg: f64, a_std: f64) -> HelicalHBond {
        let (x_from, x_to) = plus_minus_3_sigma(x_avg, x_std);
        let (r_from, r_to) = plus_minus_3_sigma(r_avg, r_std);
        let (a_from, a_to) = plus_minus_3_sigma(a_avg, a_std);

        return HelicalHBond {
            y_center, z_center,
            function: HBondFunction{x_from, x_to, r_from, r_to, a_from, a_to, smooth: 0.2},
        };
    }
    pub fn evaluate(&self,ca_prev: &Vec3, ca_i: &Vec3, ca_next: &Vec3, ca_j: &Vec3) -> f64{
        // ---------- find x,r,a coordinates of ca_j in the LCS of ca_i
        let (x, r, a) = compute_x_r_a(&ca_prev, &ca_i, &ca_next, &ca_j, self.y_center, self.z_center);
        // ---------- calculate the bond probability
        return self.function.evaluate(x, r, a);
    }
}

struct StrandHBond {
    y_center: f64,
    z_center: f64,
    function_a_lt_5: HBondFunction,
    function_a_ge_5: HBondFunction,
}

impl StrandHBond {
    pub fn new(y_center: f64, z_center: f64, x_avg_lt: f64, x_std_lt: f64, r_avg_lt: f64, r_std_lt: f64, a_avg_lt: f64, a_std_lt: f64,
               x_avg_ge: f64, x_std_ge: f64, r_avg_ge: f64, r_std_ge: f64, a_avg_ge: f64, a_std_ge: f64) -> StrandHBond {

        let (x_from, x_to) = plus_minus_3_sigma(x_avg_lt, x_std_lt);
        let (r_from, r_to) = plus_minus_3_sigma(r_avg_lt, r_std_lt);
        let (a_from, a_to) = plus_minus_3_sigma(a_avg_lt, a_std_lt);
        let function_a_lt_5 = HBondFunction { x_from, x_to, r_from, r_to, a_from, a_to, smooth: 0.2 };
        let (x_from, x_to) = plus_minus_3_sigma(x_avg_ge, x_std_ge);
        let (r_from, r_to) = plus_minus_3_sigma(r_avg_ge, r_std_ge);
        let (a_from, a_to) = plus_minus_3_sigma(a_avg_ge, a_std_ge);
        let function_a_ge_5 = HBondFunction { x_from, x_to, r_from, r_to, a_from, a_to, smooth: 0.2 };

        return StrandHBond {
            y_center, z_center,
            function_a_lt_5, function_a_ge_5
        };
    }
    pub fn evaluate(&self,ca_prev: &Vec3, ca_i: &Vec3, ca_next: &Vec3, ca_j: &Vec3) -> f64{
        // ---------- find x,r,a coordinates of ca_j in the LCS of ca_i
        let (x, r, a) = compute_x_r_a(&ca_prev, &ca_i, &ca_next, &ca_j, self.y_center, self.z_center);
        // ---------- calculate the bond probability
        if a >= 5.0 { return self.function_a_ge_5.evaluate(x, r, a); }
        else {return self.function_a_lt_5.evaluate(x, r, a);}
    }
}


fn plus_minus_3_sigma(avg: f64, std: f64) -> (f64, f64) {
    return (avg - 3.0 * std, avg + 3.0 * std);
}

enum HBondSense {
    None,
    Parallel,
    Antiparallel
}
fn sense_h_bond(ca_a_prev: &Vec3, ca_a: &Vec3, ca_a_next: &Vec3,
                ca_b_prev: &Vec3, ca_b: &Vec3, ca_b_next: &Vec3) -> HBondSense {

    let d = ca_a.distance_square_to(ca_b);
    if d < MIN_HBOND_DISTANCE_SQ || d > MAX_HBOND_DISTANCE_SQ { return  HBondSense::None; }

    let d_ap_bp = ca_a_prev.distance_square_to(ca_b_prev);
    let d_an_bn = ca_a_next.distance_square_to(ca_b_next);
    if d_ap_bp >= MIN_HBOND_DISTANCE_SQ && d_ap_bp <= MAX_HBOND_DISTANCE_SQ &&
        d_an_bn >= MIN_HBOND_DISTANCE_SQ && d_an_bn <= MAX_HBOND_DISTANCE_SQ { return  HBondSense::Parallel; }

    let d_ap_bn = ca_a_prev.distance_square_to(ca_b_next);
    let d_an_bp = ca_a_next.distance_square_to(ca_b_prev);
    if d_ap_bn >= MIN_HBOND_DISTANCE_SQ && d_ap_bn <= MAX_HBOND_DISTANCE_SQ &&
        d_an_bp >= MIN_HBOND_DISTANCE_SQ && d_an_bp <= MAX_HBOND_DISTANCE_SQ { return  HBondSense::Antiparallel; }

    return HBondSense::None;
}



/// Calculates x, r, alpha coordinates of a CA atom from its Cartesian coordinates
///
/// This function computes internal coordinates of the ``ca_j`` atom in the Local coordinate System
/// formed by the three atoms: ``ca_prev``, ``ca_i`` and ``ca_next``. Subsequently it computes
/// the `x``, ``r`` and ``alpha`` corrected coordinates of that point
fn compute_x_r_a(ca_prev: &Vec3, ca_i: &Vec3, ca_next: &Vec3, ca_j: &Vec3,
                 correct_y: f64, correct_z: f64) -> (f64, f64, f64) {

    let rt = Rototranslation::by_three_atoms(&ca_prev, &ca_i, &ca_next);
    let mut ca_j_local = rt.apply(&ca_j);
    // println!("xyz: {} {} {}", ca_j_local.x, ca_j_local.y, ca_j_local.z);

    let x = ca_j_local.x;                                        // --- X coordinate
    ca_j_local.y -= correct_y;                                        // --- correct for the center of Y-Z circle
    ca_j_local.z -= correct_z;
    // if (ca_j.z > 4.0) return false;                          // --- "stacking" strands is not allowed!
    let r = (ca_j_local.y * ca_j_local.y + ca_j_local.z * ca_j_local.z).sqrt();    // --- polar coordinates from Y-Z
    let mut a = ca_j_local.z.atan2(ca_j_local.y);                // --- angle should be [-pi,pi]
    if a < 0.0 { a += 2.0 * std::f64::consts::PI; }               // --- take care the angle is in the range [0,2pi]
    if a < std::f64::consts::PI / 2.0 { a += 2.0 * std::f64::consts::PI; }

    // println!("{x} {r} {a}");
    return (x, r, a);
}


#[cfg(test)]
mod tests {
    use bioshell_pdb::assert_delta;
    use super::*;

    #[test]
    fn test_x_r_a_computations() {
        let vip = Vec3::new(-3.0, 0.0, 0.0);
        let vi = Vec3::new(0.0, 2.3, 0.0);
        let vin = Vec3::new(3.0, 0.0, 0.0);
        let mut vjp = Vec3::new(-3.0, 0.0, 5.0);
        let mut vj = Vec3::new(0.0, 2.3, 5.0);
        let mut vjn = Vec3::new(3.0, 0.0, 5.0);

        let (x, r, a) = compute_x_r_a(&vip, &vi, &vin, &vj, 0.0, 0.0);
        assert_delta!(x, 0.0, 0.00001);
        assert_delta!(r, 5.0, 0.00001);
        assert_delta!(a, 3.14159, 0.0001);
    }
}