use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::time::Instant;
use flate2::read::GzDecoder;
use log::info;
use tar::Archive;
use reqwest::blocking::get;

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

#[derive(Debug)]
pub struct Node {
    pub tax_id: u32,
    pub parent_tax_id: u32,
    pub rank: Rank,
    pub name: String,
}

#[derive(Debug)]
pub struct Taxonomy {
    nodes: Vec<Node>,
    taxid_to_index: HashMap<u32, usize>,
    name_to_taxid: HashMap<String, u32>,
}

impl Taxonomy {

    pub fn species(&self, taxid: u32) -> Option<&Node> {
        let mut current = self.taxid_to_index.get(&taxid).cloned();
        while let Some(idx) = current {
            let node = &self.nodes[idx];
            if node.rank == Rank::Species {
                return Some(node);
            }
            current = self.taxid_to_index.get(&node.parent_tax_id).cloned();
        }
        None
    }

    pub fn lineage(&self, taxid: u32) -> Vec<&Node> {
        let mut lineage = Vec::new();
        let mut current = self.taxid_to_index.get(&taxid).cloned();
        while let Some(idx) = current {
            let node = &self.nodes[idx];
            lineage.push(node);
            if node.tax_id == node.parent_tax_id {
                break;
            }
            current = self.taxid_to_index.get(&node.parent_tax_id).cloned();
        }
        lineage.reverse();
        lineage
    }

    pub fn taxid(&self, name: &str) -> Option<u32> {
        self.name_to_taxid.get(name).cloned()
    }

    pub fn rank(&self, tax_id: u32, rank: Rank) -> Option<&Node> {
        let mut current = self.taxid_to_index.get(&tax_id).cloned();
        while let Some(idx) = current {
            let node = &self.nodes[idx];
            if node.rank == rank {
                return Some(node);
            }
            if node.tax_id == node.parent_tax_id {
                break; // reached root
            }
            current = self.taxid_to_index.get(&node.parent_tax_id).cloned();
        }
        None
    }

    pub fn load_from_tar_gz<P: AsRef<Path>>(path: P) -> Result<Self, Box<dyn Error>> {
        let start = Instant::now();

        let file = File::open(path)?;
        let gz = GzDecoder::new(file);
        let mut archive = Archive::new(gz);

        let mut nodes_raw: Vec<(u32, u32, Rank)> = Vec::new();
        let mut names_raw: HashMap<u32, String> = HashMap::new();

        for entry in archive.entries()? {
            let entry = entry?;
            let path = entry.path()?;

            if path.ends_with("nodes.dmp") {
                let reader = BufReader::new(entry);
                for line in reader.lines() {
                    let line = line?;
                    let fields: Vec<&str> = line.split("\t|\t").map(|s| s.trim_matches('|')).collect();
                    if fields.len() >= 3 {
                        let tax_id = fields[0].trim().parse()?;
                        let parent_tax_id = fields[1].trim().parse()?;
                        let rank = Rank::from_str(fields[2].trim());
                        nodes_raw.push((tax_id, parent_tax_id, rank));
                    }
                }
            } else if path.ends_with("names.dmp") {
                let reader = BufReader::new(entry);
                for line in reader.lines() {
                    let line = line?;
                    let fields: Vec<&str> = line.split("\t|\t").map(|s| s.trim_matches('|')).collect();
                    if fields.len() >= 4 && fields[3].trim() == "scientific name" {
                        let tax_id = fields[0].trim().parse()?;
                        let name = fields[1].trim().to_string();
                        names_raw.insert(tax_id, name);
                    }
                }
            }
        }

        let mut nodes = Vec::new();
        let mut taxid_to_index = HashMap::new();
        let mut name_to_taxid = HashMap::new();

        for (i, (tax_id, parent_tax_id, rank)) in nodes_raw.into_iter().enumerate() {
            let name = names_raw.get(&tax_id).cloned().unwrap_or_else(|| "".to_string());
            if !name.is_empty() {
                name_to_taxid.insert(name.clone(), tax_id);
            }
            nodes.push(Node { tax_id, parent_tax_id, rank, name });
            taxid_to_index.insert(tax_id, i);
        }

        info!("Loaded taxonomy in {:.2?}", start.elapsed());

        Ok(Taxonomy { nodes, taxid_to_index, name_to_taxid })
    }

    #[allow(dead_code)]
    fn download_taxdump_to_file(url: &str, output_path: &Path) -> Result<(), Box<dyn Error>> {
        let start = Instant::now();

        let mut file = File::create(output_path)?;
        let response = get(url)?;
        let content = response.bytes()?;
        file.write_all(&content)?;

        let elapsed = start.elapsed();
        info!("Downloaded taxonomy dump in {:.2?}", elapsed);

        Ok(())
    }
}
