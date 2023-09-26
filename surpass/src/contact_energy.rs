use crate::{MoveProposal, SurpassAlphaSystem, SurpassEnergy};



struct CaContactEnergy {
    e_cont: f64,
    e_rep: f64,
    r_max: f64,
    r_min: f64,
    r_rep: f64,
    i_rep_2: f64,
    i_min_2: f64,
    i_max_2: f64,
}

impl CaContactEnergy {
    pub fn new(system: &SurpassAlphaSystem, e_rep: f64, e_cont: f64, r_rep: f64, r_min: f64, r_max: f64) -> CaContactEnergy {
        let mut i_rep_2 = system.real_to_int(r_rep) as f64;
        i_rep_2 *= i_rep_2;
        let mut i_min_2 = system.real_to_int(r_min) as f64;
        i_min_2 *= i_min_2;
        let mut i_max_2 = system.real_to_int(r_max) as f64;
        i_max_2 *= i_max_2;

        CaContactEnergy{ e_cont, e_rep, r_max, r_min, r_rep, i_rep_2, i_min_2, i_max_2 }
    }
}

macro_rules! diff_squared_or_continue {
    ($xi: expr, $xj: expr, $sum: expr, $cutoff2: expr) => {
        let mut dx: f64 = ($xi - $xj) as f64;
        dx *= dx;
        $sum += dx;
        if dx >= $cutoff2 { continue }
    }
}

// macro_rules! energy_for_distance {
//     ($self: expr, $d2: expr, $e_tot: expr) => {
//         if d2 < $self.i_rep_2 { $e_tot+= $self.e_rep; }
//     }
// }

macro_rules! energy_for_distance {
    ($self: expr, $d2: expr) => {
        if $d2 < $self.i_rep_2 {  $self.e_rep }
        else { 0.0 }
    }
}
impl SurpassEnergy for CaContactEnergy {

    fn evaluate(&self, conf: &SurpassAlphaSystem) -> f64 {
        todo!()
    }

    fn evaluate_delta<const N: usize>(&self, conf: &SurpassAlphaSystem, move_prop: &MoveProposal<N>) -> f64 {

        let mut en_chain = 0.0;
        let mut en_proposed = 0.0;
        let mut sum_d2: f64;
        let mut i_chain = move_prop.first_moved_pos;
        for i_moved in 0..N {
            for i_partner in 0..(move_prop.first_moved_pos-1) {
                en_proposed += 'energy_after: {
                    let dx = (move_prop.cax[i_moved] - conf.cax[i_partner]) as f64;
                    sum_d2 = dx * dx;
                    if sum_d2 > self.i_max_2 { break 'energy_after 0.0 }
                    let dy = (move_prop.cay[i_moved] - conf.cay[i_partner]) as f64;
                    sum_d2 += dy * dy;
                    if sum_d2 > self.i_max_2 { break 'energy_after 0.0 }
                    let dz = (move_prop.caz[i_moved] - conf.caz[i_partner]) as f64;
                    sum_d2 += dz * dz;
                    energy_for_distance!(self, sum_d2)
                };
                en_chain += 'energy_before: {
                    let dx = (conf.cax[i_chain] - conf.cax[i_partner]) as f64;
                    sum_d2 = dx * dx;
                    if sum_d2 > self.i_max_2 { break 'energy_before 0.0 }
                    let dy = (conf.cay[i_chain] - conf.cay[i_partner]) as f64;
                    sum_d2 += dy * dy;
                    if sum_d2 > self.i_max_2 { break 'energy_before 0.0 }
                    let dz = (conf.caz[i_chain] - conf.caz[i_partner]) as f64;
                    sum_d2 += dz * dz;
                    energy_for_distance!(self, sum_d2)
                };
                i_chain += 1;
            }
            // for i_partner in (move_prop.first_moved_pos+N+1)..conf.count_atoms() {
            //     sum_d2 = 0.0;
            //     diff_squared_or_continue!(move_prop.moved_cax[i_moved], conf.cax[i_partner], sum_d2, self.r_max_2);
            //     diff_squared_or_continue!(move_prop.moved_cay[i_moved], conf.cay[i_partner], sum_d2, self.r_max_2);
            //     diff_squared_or_continue!(move_prop.moved_caz[i_moved], conf.caz[i_partner], sum_d2, self.r_max_2);
            //     energy_for_distance!(self, sum_d2, en);
            // }
            // --- apply energy
        }

        return en_proposed - en_chain;
    }
}