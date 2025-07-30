use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::time::Instant;
use flate2::read::GzDecoder;
use log::{debug, info};
use tar::Archive;
use reqwest::blocking::get;
use serde::Serialize;
use crate::Rank;

/// A node of a taxonomy tree.
///
/// A node can represent a species, genus, family, order, etc.
///
/// # Examples
/// ```
/// use bioshell_taxonomy::{Taxonomy, Rank};
/// # use std::error::Error;
/// # fn main() -> Result<(), Box<dyn Error>> {
/// # let path = "./tests/test_files/test_taxdump.tar.gz";
/// let taxonomy = Taxonomy::load_from_tar_gz(&path)?;
/// let human_node = taxonomy.node(9606).unwrap();
/// assert_eq!(human_node.name, "Homo sapiens");
/// assert_eq!(human_node.rank, Rank::Species);
/// let human_node = taxonomy.node(9443).unwrap();
/// assert_eq!(human_node.name, "Primates");
/// assert_eq!(human_node.rank, Rank::Order);
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Serialize)]
pub struct Node {
    /// Taxonomy-wide unique identifier
    pub tax_id: u32,
    /// ``taxid`` identifier of the parent node
    pub parent_tax_id: u32,
    /// Rank of the node, e.g. [``Rank::Genus``](Rank::Genus) or [``Rank::Phylum``](Rank::Phylum)
    pub rank: Rank,
    /// Scientific name, e.g. "Homo sapiens"
    pub name: String,
}

/// Stores a taxonomy tree, build from a NCBI taxonomy `taxdump.tar.gz` dump.
///
/// The taxonomy data file can be downloaded from NCBI:
/// <https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz>
///
/// Having the file on your computer, you can load it into a `Taxonomy` object:
/// ```no_run
/// use bioshell_taxonomy::Taxonomy;
/// let taxonomy = Taxonomy::load_from_tar_gz("taxdump.tar.gz")
///             .expect("Can't load taxdump.tar.gz file!");
/// ```
/// which then you can use to find information about any taxid:
/// ```
/// # use bioshell_taxonomy::{Rank, Taxonomy};
/// # let taxonomy = Taxonomy::load_from_tar_gz("./tests/test_files/test_taxdump.tar.gz").expect("Can't load taxdump.tar.gz file!");
/// let human = taxonomy.node(9606).unwrap(); // returns the node for human
/// assert_eq!(human.name, "Homo sapiens");
/// assert_eq!(human.rank, Rank::Species);
/// ```
/// It's also possible to find taxonomy information by name:
/// ```
/// # use bioshell_taxonomy::{Rank, Taxonomy};
/// # let taxonomy = Taxonomy::load_from_tar_gz("./tests/test_files/test_taxdump.tar.gz").expect("Can't load taxdump.tar.gz file!");
/// if let Some(human_taxid) = taxonomy.taxid("Homo sapiens") {
///     let human_node = taxonomy.node(human_taxid).unwrap();
///     assert_eq!(human_node.rank, Rank::Species);
/// }
/// ```
///

#[derive(Debug)]
pub struct Taxonomy {
    nodes: Vec<Node>,
    taxid_to_index: HashMap<u32, usize>,
    // maps any name to its taxid, e.g. "human" -> 9606, "Homo sapiens" -> 9606
    pub(crate) name_to_taxid: HashMap<String, u32>,
}

impl Taxonomy {

    /// Returns the node for a given ``taxid``.
    ///
    /// # Examples
    /// ```
    /// # use bioshell_taxonomy::Taxonomy;
    /// # use std::error::Error;
    /// # fn main() -> Result<(), Box<dyn Error>> {
    /// # use bioshell_taxonomy::Rank;
    /// # let path = "./tests/test_files/test_taxdump.tar.gz";
    /// let taxonomy = Taxonomy::load_from_tar_gz(&path)?;
    /// let human_node = taxonomy.node(9606).unwrap();
    /// assert_eq!(human_node.name, "Homo sapiens");
    /// assert_eq!(human_node.rank, Rank::Species);
    /// # Ok(())
    /// # }
    /// ```
    pub fn node(&self, taxid: u32) -> Option<&Node> {
        let current = self.taxid_to_index.get(&taxid).cloned();
        if let Some(idx) = current {
            return Some(&self.nodes[idx]);
        } else { None }
    }

    /// Provides an iterator over all nodes in this taxonomy.
    ///
    /// The nodes provided by the iterator can be filtered, e.g by rank:
    ///
    /// # Example
    /// ```no_run
    /// use bioshell_taxonomy::{Rank, Taxonomy};
    /// let taxonomy = Taxonomy::load_from_tar_gz("./tests/test_files/test_taxdump.tar.gz").unwrap();
    /// for phylum  in taxonomy.nodes().filter(|node| node.rank == Rank::Phylum) {
    ///     println!("{}", phylum.name);
    /// }
    /// ```
    pub fn nodes(&self) -> impl Iterator<Item = &Node> {
        self.nodes.iter()
    }

    /// Iterates over all the names associated with a given `taxid`
    ///
    /// The iterator provides the scientific name, synonyms as well as common names,
    /// as they are defined in the NCBI taxonomy.
    pub fn names(&self, taxid: u32) -> impl Iterator<Item = &String> {
        self.name_to_taxid
            .iter()
            .filter(move |&(_, &id)| id == taxid)
            .map(|(name, _)| name)
    }

    /// Provides the lineage of a given taxonomy node, e.g. for a species.
    ///
    /// The lineage is returned as a vector of [`Node`]s, starting from the root to the given taxid.
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

    /// Finds a species given its name and returns `taxid`
    /// ```
    /// # use bioshell_taxonomy::{Rank, Taxonomy};
    /// # let taxonomy = Taxonomy::load_from_tar_gz("./tests/test_files/test_taxdump.tar.gz").expect("Can't load taxdump.tar.gz file!");
    /// if let Some(human_taxid) = taxonomy.taxid("Homo sapiens") {
    ///     let human_node = taxonomy.node(human_taxid).unwrap();
    ///     assert_eq!(human_node.rank, Rank::Species);
    /// }
    /// ```
    pub fn taxid(&self, name: &str) -> Option<u32> {
        debug!("Resolving name: {}", name);
        self.name_to_taxid.get(name).cloned()
    }

    /// Look for a specific rank in the lineage of a given node.
    ///
    /// This method traverses the taxonomy tree staring from the node given by taxid to the top
    /// of the tree and returns the first node that matches the given rank.
    ///
    ///
    /// # Examples
    /// ```
    /// use bioshell_taxonomy::{Taxonomy, Rank};
    /// # use std::error::Error;
    /// # fn main() -> Result<(), Box<dyn Error>> {
    /// # let path = "./tests/test_files/test_taxdump.tar.gz";
    /// let taxonomy = Taxonomy::load_from_tar_gz(&path)?;
    /// let order = taxonomy.rank(9606, Rank::Order).unwrap();
    /// assert_eq!(order.name, "Primates");
    /// assert_eq!(order.rank, Rank::Order);
    /// # Ok(())
    /// # }
    /// ```
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

    /// Create the taxonomy from a given NCBI taxonomy dump.
    ///
    /// The taxonomy data should be downloaded from NCBI website. Don't unpack the file; this method loads the whole archive.
    /// # Example
    /// ``` no_run
    /// use bioshell_taxonomy::{Taxonomy};
    /// let taxonomy = Taxonomy::load_from_tar_gz("./taxdump.tar.gz").expect("Can't load taxdump.tar.gz file!");
    /// ```
    pub fn load_from_tar_gz<P: AsRef<Path>>(path: P) -> Result<Self, Box<dyn Error>> {
        let start = Instant::now();

        let file = File::open(path)?;
        let gz = GzDecoder::new(file);
        let mut archive = Archive::new(gz);

        let mut nodes_raw: Vec<(u32, u32, Rank)> = Vec::new();
        let mut names_raw: HashMap<u32, String> = HashMap::new();
        let mut name_to_taxid = HashMap::new();
        let mut primary_names: HashMap<u32, String> = HashMap::new();

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
                    if fields.len() >= 4 {
                        let tax_id = fields[0].trim().parse()?;
                        let name = fields[1].trim().to_string();
                        let name_class = fields[3].trim();
                        // ---------- extract scientific names for Node
                        if name_class == "scientific name" {
                            primary_names.insert(tax_id, name.clone());
                        }
                        // --------- Also allow searching by these names
                        if matches!(name_class, "synonym" | "common name" | "genbank common name") {
                            names_raw.insert(tax_id, name.clone());
                            name_to_taxid.insert(name, tax_id);
                        }
                    }
                }
            }
        }

        let mut nodes = Vec::new();
        let mut taxid_to_index = HashMap::new();

        for (i, (tax_id, parent_tax_id, rank)) in nodes_raw.into_iter().enumerate() {
            let name = primary_names.get(&tax_id).cloned().unwrap_or_else(|| "".to_string());
            if !name.is_empty() {
                name_to_taxid.insert(name.clone(), tax_id);
            }

            nodes.push(Node { tax_id, parent_tax_id, rank, name });
            taxid_to_index.insert(tax_id, i);
        }

        info!("Loaded taxonomy in {:.2?}", start.elapsed());

        Ok(Taxonomy { nodes, taxid_to_index, name_to_taxid })
    }

    const TEST_TAXDUMP_URL: &'static str = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz";

    /// Download the latest NCBI taxonomy dump from the NCBI FTP server.
    ///
    /// This method downloads the [`"https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"`](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz) file
    /// to the current folder.
    pub fn download_from_ncbi<P: AsRef<Path>>(output_path: P) -> Result<(), Box<dyn Error>> {
        let start = Instant::now();

        let mut file = File::create(output_path)?;
        let response = get(Self::TEST_TAXDUMP_URL)?;
        let content = response.bytes()?;
        file.write_all(&content)?;

        let elapsed = start.elapsed();
        info!("Downloaded taxonomy dump in {:.2?}", elapsed);

        Ok(())
    }
}
