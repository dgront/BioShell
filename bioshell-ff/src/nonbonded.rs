use std::ops::Range;

use crate::ff::Energy;
use crate::Coordinates;

pub struct SimpleContact {
    /// repulsion distance
    r_c: f32,
    /// energy well starts
    r_from: f32,
    /// energy well end
    r_to: f32,
    /// energy value
    en: f32,
    /// repulsion value
    rep: f32,

    sequence_separation: usize,

    r_c_sq: f32,
    r_from_sq: f32,
    r_to_sq: f32,
}

macro_rules! pairwise_contact_kernel {
    ($x:expr,$y:expr,$z:expr,$chain:expr,$i:expr,$r_to_2:expr,$r_from_2:expr,$r_rep_2:expr,$e_c:expr,$e_rep:expr,$en:expr)=>{

        let mut d = $chain[$i].x - $x;
        let mut d2 = d*d;
        if d2 > $r_to_2 { continue; }
        d = $chain[$i].y - $y;
        d2 += d*d;
        if d2 > $r_to_2 { continue; }
        d = $chain[$i].z - $z;
        d2 += d*d;
        if d2 > $r_to_2 { continue; }
        if d2 < $r_rep_2 { $en += $e_rep as f64}
        if d2 > $r_from_2 { $en += $e_c as f64 }
    }
}

macro_rules! pairwise_contact_kernel2 {
    ($x:expr, $y:expr, $z:expr, $chain:expr, $i:expr, $self:expr, $en:expr) => {

        let mut d = $chain[$i].x - $x;
        let mut d2 = d*d;
        if d2 > $self.r_to_sq { continue; }
        d = $chain[$i].y - $y;
        d2 += d*d;
        if d2 > $self.r_to_sq { continue; }
        d = $chain[$i].z - $z;
        d2 += d*d;
        if d2 > $self.r_to_sq { continue; }
        if d2 < $self.r_c_sq { $en += $self.rep as f64}
        if d2 > $self.r_from_sq { $en += $self.en as f64 }
    }
}

impl SimpleContact {
    pub fn new(r_c: f32, r_from: f32, r_to: f32, rep: f32, en: f32, separation: usize) -> SimpleContact {
        SimpleContact { r_c, r_from, r_to, rep, en, sequence_separation: separation,
            r_c_sq: r_c * r_c, r_from_sq: r_from * r_from, r_to_sq: r_to * r_to
        }
    }

    pub fn pair_energy(&self, system: &Coordinates, ipos: usize, jpos: usize) -> f32 {
        if ipos > jpos && ipos - jpos <= self.sequence_separation { return 0.0; }
        if ipos < jpos && jpos - ipos <= self.sequence_separation { return 0.0; }
        let d = system.distance_square(ipos, jpos).sqrt();
        if d > self.r_to { return 0.0; }
        if d < self.r_c { return self.rep }
        if d > self.r_from { return self.en; } else { return 0.0 }
    }

    fn each_vs_each_energy(&self, system: &Coordinates, moved: &Range<usize>) ->f64 {

        let mut en: f64 = 0.0;
        for ipos in moved.start..moved.end {
            let xi: f32 = system[ipos].x;
            let yi: f32 = system[ipos].y;
            let zi: f32 = system[ipos].z;

            for jpos in (ipos+1+self.sequence_separation)..moved.end + 1 {
                pairwise_contact_kernel!(xi, yi, zi, system, jpos, self.r_to_sq, self.r_from_sq, self.r_c_sq, self.en, self.rep, en);
            }
        }
        return en;
    }
}

impl Energy for SimpleContact {

    fn energy(&self, system: &Coordinates) -> f64 {
        let mut en:f64 = 0.0;
        for i in 0..system.size() {
            en += self.energy_by_pos(system, i);
        }

        return en / 2.0;
    }

    fn energy_by_pos(&self, chain: &Coordinates, pos:usize) -> f64 {
        let mut en: f64 = 0.0;

        let x: f32 = chain[pos].x;
        let y: f32 = chain[pos].y;
        let z: f32 = chain[pos].z;

        if pos > self.sequence_separation {
            for i in 0..(pos - self.sequence_separation) {
                // pairwise_contact_kernel!(x, y, z, chain, i, self.r_to_sq, self.r_from_sq, self.r_c_sq, self.en, self.rep, en);
                pairwise_contact_kernel2!(x, y, z, chain, i, self, en);
            }
        }
        // for the chain of 10 atoms and sep=1, we score when pos is at most 7 (7 vs 9)
        if pos < chain.size() - 1 - self.sequence_separation {
            for i in (pos + 1 + self.sequence_separation)..chain.size() {
                pairwise_contact_kernel!(x, y, z, chain, i, self.r_to_sq, self.r_from_sq, self.r_c_sq, self.en, self.rep, en);
            }
        }

        return en;
    }

    fn delta_energy_by_range(&self, old: &Coordinates, moved: &Range<usize>, new: &Coordinates) -> f64 {

        let mut en_old: f64 = 0.0;
        let mut en_new: f64 = 0.0;

        for jpos in moved.start..moved.end + 1 {
            // --- x, y, z of a moved atom are different between new and old systems
            let xo: f32 = old[jpos].x;
            let yo: f32 = old[jpos].y;
            let zo: f32 = old[jpos].z;
            let xn: f32 = new[jpos].x;
            let yn: f32 = new[jpos].y;
            let zn: f32 = new[jpos].z;

            // --- Energy upstream the moved range; e.g. for sep=1 we calculate energy here if moved 2 or greater
            if moved.start > self.sequence_separation {
                for ipos in 0..(moved.start - self.sequence_separation) {
                    pairwise_contact_kernel!(xo, yo, zo, old, ipos, self.r_to_sq, self.r_from_sq, self.r_c_sq, self.en, self.rep, en_old);
                }
                for ipos in 0..(moved.start - self.sequence_separation) {
                    pairwise_contact_kernel!(xn, yn, zn, new, ipos, self.r_to_sq, self.r_from_sq, self.r_c_sq, self.en, self.rep, en_new);
                }
            }
            // --- Energy downstream the moved range, e.g. for N=10 and sep=1 we calculate here if moved at most 7 (7 vs 9)
            if moved.end < old.size() - self.sequence_separation -1 {
                let start = moved.end + 1 + self.sequence_separation;
                for ipos in start..old.size() {
                    pairwise_contact_kernel!(xo, yo, zo, old, ipos, self.r_to_sq, self.r_from_sq, self.r_c_sq, self.en, self.rep, en_old);
                }
                for ipos in start..old.size() {
                    pairwise_contact_kernel!(xn, yn, zn, new, ipos, self.r_to_sq, self.r_from_sq, self.r_c_sq, self.en, self.rep, en_new);
                }
            }
        }
        // --- Energy within the moved range
        en_old += self.each_vs_each_energy(old,moved);
        en_new += self.each_vs_each_energy(new,moved);

        return en_new - en_old;
    }
}