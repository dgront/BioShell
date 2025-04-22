#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum Rank {
    Unclassified = 0,
    Species = 1,
    Genus = 2,
    Family = 3,
    Order = 4,
    Class = 5,
    Phylum = 6,
    Kingdom = 7,
    Superkingdom = 8,
    Other = 255,
}

impl Rank {
    pub fn from_str(s: &str) -> Self {
        let s_lower = s.to_lowercase();
        match s_lower.as_str() {
            "species" => Rank::Species,
            "genus" => Rank::Genus,
            "family" => Rank::Family,
            "order" => Rank::Order,
            "class" => Rank::Class,
            "phylum" => Rank::Phylum,
            "kingdom" => Rank::Kingdom,
            "superkingdom" => Rank::Superkingdom,
            "domain" => Rank::Superkingdom,
            _ => Rank::Other,
        }
    }
}

use std::cmp::Ordering;

impl PartialOrd for Rank {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Rank {
    fn cmp(&self, other: &Self) -> Ordering {
        (*self as u8).cmp(&(*other as u8))
    }
}