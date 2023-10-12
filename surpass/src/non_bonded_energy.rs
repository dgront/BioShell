use crate::{MoveProposal, NonBondedEnergyKernel, SurpassAlphaSystem, SurpassEnergy};


pub struct NonBondedEnergy<E: NonBondedEnergyKernel> {
    i_cutoff_2: f64,
    energy_kernel: E
}

impl<E: NonBondedEnergyKernel> NonBondedEnergy<E> {
    pub fn new(system: &SurpassAlphaSystem, energy_kernel: E) -> NonBondedEnergy<E> {
        let mut i_rep_2 = system.real_to_int(energy_kernel.distance_cutoff()) as f64;
        i_rep_2 *= i_rep_2;

        NonBondedEnergy{ i_cutoff_2: i_rep_2, energy_kernel}
    }
}

macro_rules! pairwise_energy {
    ($self: expr, $i: expr, $i_coords: expr, $j: expr, $j_coords: expr, $en_total: expr) => {
        $en_total += 'energy_after: {
            let dx = ($i_coords.cax[$i].wrapping_sub($j_coords.cax[$j])) as f64;
            let mut sum_d2 = dx * dx;
            if sum_d2 > $self.i_cutoff_2 { break 'energy_after 0.0 }
            let dy = ($i_coords.cay[$i].wrapping_sub($j_coords.cay[$j])) as f64;
            sum_d2 += dy * dy;
            if sum_d2 > $self.i_cutoff_2 { break 'energy_after 0.0 }
            let dz = ($i_coords.caz[$i].wrapping_sub($j_coords.caz[$j])) as f64;
            sum_d2 += dz * dz;
            $self.energy_kernel.energy_for_distance_squared(sum_d2)
        };
    }
}

impl<E: NonBondedEnergyKernel> SurpassEnergy for NonBondedEnergy<E> {
    fn evaluate(&self, conf: &SurpassAlphaSystem) -> f64 {
        let mut e_total = 0.0;
        for i in 1..conf.count_atoms() as i32 {
            for j in 0..i {
                pairwise_energy!(self, i as usize, conf, j as usize, conf, e_total);
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
                pairwise_energy!(self, i_partner as usize, conf, i_moved, move_prop, en_proposed);
                pairwise_energy!(self, i_partner as usize, conf, i_chain as usize, conf, en_chain);
                // eprintln!("{} {} {}   {} {}", i_partner, i_moved, i_chain, en_chain, en_proposed);
            }
            i_chain += 1;
        }
        let mut i_chain = move_prop.first_moved_pos as i32;
        for i_moved in 0..N {
            for i_partner in (move_prop.first_moved_pos + N) as i32 ..conf.count_atoms() as i32 {
                pairwise_energy!(self, i_partner as usize, conf, i_moved, move_prop, en_proposed);
                pairwise_energy!(self, i_partner as usize, conf, i_chain as usize, conf, en_chain);
                // eprintln!("{} {} {}   {} {}", i_partner, i_moved, i_chain, en_chain, en_proposed);
            }
            i_chain += 1;
        }

        return en_proposed - en_chain;
    }
}