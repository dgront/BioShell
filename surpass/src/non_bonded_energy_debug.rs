use crate::{MoveProposal, NonBondedEnergyKernel, SurpassAlphaSystem, SurpassEnergy};
use crate::non_bonded_energy::pairwise_energy;
use crate::SurpassAtomTypes;


pub struct NonBondedEnergyDebug<E: NonBondedEnergyKernel> {
    i_cutoff_2: f64,
    energy_kernel: E,
    global_energy_map: Vec<Vec<f64>>,
    delta_energy_map: Vec<Vec<f64>>,
}

impl<E: NonBondedEnergyKernel> NonBondedEnergyDebug<E> {
    pub fn new(system: &SurpassAlphaSystem, energy_kernel: E) -> NonBondedEnergyDebug<E> {
        let mut i_rep_2 = system.real_to_int(energy_kernel.distance_cutoff()) as f64;
        i_rep_2 *= i_rep_2;
        let n_atoms = system.count_residues();
        let global_energy_map = vec![vec![0.0; n_atoms]; n_atoms];
        let delta_energy_map = vec![vec![0.0; n_atoms]; n_atoms];
        NonBondedEnergyDebug{ i_cutoff_2: i_rep_2, energy_kernel, global_energy_map, delta_energy_map }
    }

    pub fn evaluate(&mut self, conf: &SurpassAlphaSystem) -> f64 {

        for i in 0..conf.count_residues() {
            self.global_energy_map[i].fill(0.0);
        }

        let mut e_total = 0.0;
        for i in 1..conf.count_residues() as i32 {
            for j in 0..i {
                let old_en = e_total;
                pairwise_energy!(self, i as usize, conf, j as usize, conf, e_total);
                if e_total != old_en {
                    self.global_energy_map[j as usize][i as usize] = old_en;
                    self.global_energy_map[i as usize][j as usize] = e_total;
                }
            }
        }
        return e_total;
    }

    pub fn evaluate_delta<const N_RESIDUES: usize, const N_ATOMS: usize>(&mut self, conf: &SurpassAlphaSystem, move_prop: &MoveProposal<N_RESIDUES, N_ATOMS>) -> f64 {

        for i in 0..conf.count_residues() {
            self.delta_energy_map[i].fill(0.0);
        }

        let mut en_chain = 0.0;
        let mut en_proposed = 0.0;
        let mut i_chain = move_prop.first_moved_pos as i32;
        for i_moved in 0..N_RESIDUES {
            for i_partner in 0..move_prop.first_moved_pos as i32 {
                let mut ep = 0.0;
                pairwise_energy!(self, i_partner as usize, conf, i_moved, move_prop, ep);
                en_proposed += ep;
                let mut ec = 0.0;
                pairwise_energy!(self, i_partner as usize, conf, i_chain as usize, conf, ec);
                en_chain += ec;
                if ep != ec {
                    let i = i_partner.max(i_moved as i32);
                    let j = i_partner.min(i_moved as i32);
                    self.global_energy_map[j as usize][i as usize] = ec;
                    self.global_energy_map[i as usize][j as usize] = ep;
                }
            }
            i_chain += 1;
        }
        let mut i_chain = move_prop.first_moved_pos as i32;
        for i_moved in 0..N_RESIDUES {
            for i_partner in (move_prop.first_moved_pos + N_RESIDUES) as i32 ..conf.count_residues() as i32 {
                let mut ep = 0.0;
                pairwise_energy!(self, i_partner as usize, conf, i_moved, move_prop, ep);
                en_proposed += ep;
                let mut ec = 0.0;
                pairwise_energy!(self, i_partner as usize, conf, i_chain as usize, conf, ec);
                en_chain += ec;
                if ep != ec {
                    let i = i_partner.max(i_moved as i32);
                    let j = i_partner.min(i_moved as i32);
                    self.global_energy_map[j as usize][i as usize] = ec;
                    self.global_energy_map[i as usize][j as usize] = ep;
                }
            }
            i_chain += 1;
        }

        return en_proposed - en_chain;
    }

    fn show_matrix(m: &Vec<Vec<f64>>) {
        let n = m.len();
        eprint!("    ");
        for i in 0..n { eprint!(" {:3}",i); }
        eprintln!();
        for i in 0..n {
            eprint!(" {:3}",i);
            for j in 0..n {
                eprint!(" {:3}", m[i][j]);
            }
            eprintln!();
        }
        eprintln!();
    }

    pub fn report<const N_RESIDUES: usize, const N_ATOMS: usize>(&mut self, conf: &mut SurpassAlphaSystem, move_prop: &MoveProposal<N_RESIDUES, N_ATOMS>) {
        eprintln!("Moved range: [{},{}]", move_prop.first_moved_pos, move_prop.first_moved_pos + N_RESIDUES - 1);
        eprintln!("Global energy:");
        Self::show_matrix(&self.global_energy_map);
        eprintln!("Local energy:");
        Self::show_matrix(&self.delta_energy_map);
        let mut diff: Vec<(usize, usize, f64, f64)> = vec!();
        for i in 0..self.delta_energy_map.len() {
            for j in 0..self.delta_energy_map.len() {
                let dglobal = self.global_energy_map[i][j] - self.global_energy_map[j][i];
                let dlocal = self.delta_energy_map[i][j] - self.delta_energy_map[j][i];
                if (dglobal - dlocal).abs() > 0.001 { diff.push((i, j, conf.distance(i,j), 0.0)) }
            }
        }
        conf.to_pdb_file("before.pdb", false);
        move_prop.apply(conf);
        eprintln!("Differences (i, j, dist_before, dist_after):");
        for dif in &mut diff {
            let (ii, jj, d_before, mut d_after) = dif;
            d_after= conf.distance(*ii, *jj);
            eprintln!("{} {} {} {}", ii, jj, d_before, d_after);
        }
        conf.to_pdb_file("after.pdb", false);
    }
}