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

    r_c_sq: f32,
    r_from_sq: f32,
    r_to_sq: f32,

}

impl SimpleContact {
    pub fn new(r_c: f32, r_from: f32, r_to: f32, rep: f32, en: f32) -> SimpleContact {
        SimpleContact { r_c, r_from, r_to, rep, en, r_c_sq: r_c * r_c, r_from_sq: r_from * r_from, r_to_sq: r_to * r_to }
    }
}


macro_rules! pairwise_contact_kernel {
    ($x:expr,$y:expr,$z:expr,$chain:expr,$i:expr,$A2:expr,$B2:expr,$C2:expr,$e_c:expr,$e_rep:expr,$en:expr)=>{
        let mut d = $chain[$i].x - $x;
        let mut d2 = d*d;
        if d2 > $A2 { continue; }
        d = $chain[$i].y - $y;
        d2 += d*d;
        if d2 > $A2 { continue; }
        d = $chain[$i].z - $z;
        d2 += d*d;
        if d2 < $C2 { $en += $e_rep  as f64; }
        else {
            if d2 > $B2 {
                $en += $e_c as f64;
            }
        }
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

        if pos > 1 {
            for i in 0..pos - 1 {
                pairwise_contact_kernel!(x, y, z, chain, i, self.r_to_sq, self.r_from_sq, self.r_c_sq,
                    self.en, self.rep, en);
            }
        }
        if pos < chain.size() - 2 {
            for i in pos + 2..chain.size() {
                pairwise_contact_kernel!(x, y, z, chain, i, self.r_to_sq, self.r_from_sq, self.r_c_sq,
                    self.en, self.rep, en);
            }
        }
        return en;
    }

    fn delta_energy_by_pos(&self, old: &Coordinates, pos: usize, new: &Coordinates) -> f64 {

        return 0.0;
    }
}