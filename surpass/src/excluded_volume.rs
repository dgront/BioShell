use crate::{MoveProposal, SurpassAlphaSystem, SurpassEnergy};


pub struct ExcludedVolume {
    e_penalty: f64,
    r_rep: f64,
    i_rep_2: f64,
}

impl ExcludedVolume {
    pub fn new(system: &SurpassAlphaSystem, r_rep: f64, e_rep: f64) -> ExcludedVolume {
        let mut i_rep_2 = system.real_to_int(r_rep) as f64;
        i_rep_2 *= i_rep_2;

        ExcludedVolume{ e_penalty: e_rep, r_rep, i_rep_2}
    }

    /// Returns the excluded volume distance that is currently used
    pub fn repulsion_cutoff(&self) -> f64 { self.r_rep }

    /// Returns the excluded volume penalty value that is currently used
    pub fn repulsion_energy(&self) -> f64 { self.e_penalty }

    /// Sets the new value for the excluded volume distance.
    ///
    /// A reference to the simulated `system` must be provided since the [`ExcludedVolume`] class
    /// needs to convert a double distance value to its integer representation.
    pub fn set_repulsion_cutoff(&mut self, system: &SurpassAlphaSystem, rcut: f64) {
        self.r_rep = rcut;
        self.i_rep_2 = system.real_to_int(rcut) as f64 * system.real_to_int(rcut) as f64;;
    }
}

macro_rules! excluded_volume_kernel {
    ($self: expr, $i: expr, $i_coords: expr, $j: expr, $j_coords: expr, $en_total: expr) => {
        $en_total += 'energy_after: {
            let dx = ($i_coords.cax[$i].wrapping_sub($j_coords.cax[$j])) as f64;
            let mut sum_d2 = dx * dx;
            if sum_d2 > $self.i_rep_2 { break 'energy_after 0.0 }
            let dy = ($i_coords.cay[$i].wrapping_sub($j_coords.cay[$j])) as f64;
            sum_d2 += dy * dy;
            if sum_d2 > $self.i_rep_2 { break 'energy_after 0.0 }
            let dz = ($i_coords.caz[$i].wrapping_sub($j_coords.caz[$j])) as f64;
            sum_d2 += dz * dz;

            if sum_d2 < $self.i_rep_2 { $self.e_penalty }
            else { 0.0 }
        };
    }
}

impl SurpassEnergy for ExcludedVolume {
    fn evaluate(&self, conf: &SurpassAlphaSystem) -> f64 {
        let mut e_total = 0.0;
        for i in 1..conf.count_atoms() as i32 {
            for j in 0..i {
                excluded_volume_kernel!(self, i as usize, conf, j as usize, conf, e_total);
                // println!("{} {} {}", i, j, e_total);
            }
        }
        return e_total;
    }

    fn evaluate_delta<const N: usize>(&self, conf: &SurpassAlphaSystem, move_prop: &MoveProposal<N>) -> f64 {

        let mut en_chain = 0.0;
        let mut en_proposed = 0.0;
        let mut i_chain = move_prop.first_moved_pos as i32;
        for i_moved in 0..N {
            for i_partner in 0..move_prop.first_moved_pos as i32 {
                excluded_volume_kernel!(self, i_partner as usize, conf, i_moved, move_prop, en_proposed);
                excluded_volume_kernel!(self, i_partner as usize, conf, i_chain as usize, conf, en_chain);
                // println!("{} {} {}   {} {}", i_partner, i_moved, i_chain, en_chain, en_proposed);
            }
            i_chain += 1;
        }
        let mut i_chain = move_prop.first_moved_pos as i32;
        for i_moved in 0..N {
            for i_partner in (move_prop.first_moved_pos + N) as i32 ..conf.count_atoms() as i32 {
                excluded_volume_kernel!(self, i_partner as usize, conf, i_moved, move_prop, en_proposed);
                excluded_volume_kernel!(self, i_partner as usize, conf, i_chain as usize, conf, en_chain);
                // println!("{} {} {}   {} {}", i_partner, i_moved, i_chain, en_chain, en_proposed);
            }
            i_chain += 1;
        }

        return en_proposed - en_chain;
    }
}