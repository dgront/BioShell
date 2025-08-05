use std::error::Error;
use std::sync::Arc;
use std::time::Instant;
use aho_corasick::{AhoCorasick};
use log::{debug, info};
use crate::{Node, Taxonomy};

/// Attempts to find a species name in a given string, e.g. a sequence description
///
/// Taxonomy matcher uses Aho-Corasick to perform the search as efficiently as possible.
///
/// # Examples
/// ```
/// use bioshell_taxonomy::Taxonomy;
/// # use std::error::Error;
/// # fn main() -> Result<(), Box<dyn Error>> {
/// # use bioshell_taxonomy::{Rank, TaxonomyMatcher};
/// let path = "./tests/test_files/test_taxdump.tar.gz";
/// let taxonomy = Taxonomy::load_from_tar_gz(&path)?;
/// let matcher = TaxonomyMatcher::from_taxonomy(taxonomy)?;
///
/// let is_it_mouse = matcher.find("Mus musculus");
/// assert!(is_it_mouse.is_some());
/// assert_eq!(is_it_mouse.unwrap(), 10090);                // 10090 is the taxid for mouse
///
/// let is_it_mouse = matcher.find("house mouse");          // works also with common names, if listed by NCBI
/// let borrowed_taxonomy = matcher.taxonomy();             // The original taxonomy has been consumed, but we can borrow it again
/// let mouse_node = borrowed_taxonomy.node(is_it_mouse.unwrap());
/// assert_eq!(mouse_node.unwrap().rank, Rank::Species);
/// # Ok(())
/// # }
/// ```
pub struct TaxonomyMatcher {
    taxonomy: Arc<Taxonomy>,
    matcher: AhoCorasick,
    patterns: Vec<String>, // store names for mapping back
}


impl TaxonomyMatcher {
    /// Creates a new TaxonomyMatcher object based on a given taxonomy
    pub fn from_taxonomy(taxonomy: Taxonomy)-> Result<Self, Box<dyn Error>> {

        let start = Instant::now();
        let patterns: Vec<String> = taxonomy.name_to_taxid.keys()
            .cloned().map(|s|  format!(" {} ", s))
            .collect();
        let matcher = AhoCorasick::new(&patterns)?;
        info!("Taxonomy matcher constructed in {:.2?}", start.elapsed());
        let arc = Arc::new(taxonomy);
        Ok(TaxonomyMatcher { taxonomy: Arc::clone(&arc), matcher, patterns })
    }

    /// Finds a taxid by detecting a species name in a given string.
    ///
    /// ```
    /// # use std::error::Error;
    /// # fn main() -> Result<(), Box<dyn Error>> {
    /// use bioshell_taxonomy::{Taxonomy, TaxonomyMatcher};
    /// # let path = "./tests/test_files/test_taxdump.tar.gz";
    /// let taxonomy = Taxonomy::load_from_tar_gz(&path)
    ///             .expect("Failed to load taxonomy");
    /// let matcher = TaxonomyMatcher::from_taxonomy(taxonomy)?;
    /// let is_it_mouse = matcher.find("house mouse");
    /// assert!(is_it_mouse.is_some());
    /// # Ok(())
    /// # }
    /// ```
    pub fn find(&self, description: &str) -> Option<u32> {

        let desc = Self::replace_punctuation_with_spaces(description);
        let mut matches = self.matcher.find_overlapping_iter(&desc);
        let mut matched_nodes: Vec<&Node> = vec![];
        while let Some(m) = matches.next() {
            let idx = m.pattern();
            debug!("Matching: {:?} {:?}", idx, self.patterns[idx]);
            let detected_name = &self.patterns[idx];
            let taxid = self.taxonomy.taxid(detected_name.trim()).unwrap();
            matched_nodes.push(self.taxonomy.node(taxid).unwrap());
        }
        if matched_nodes.len() > 0 {
            matched_nodes.sort_by_key(|node| node.rank);
            let mut same_rank: Vec<_> = matched_nodes.iter()
                .filter(|node| node.rank == matched_nodes[0].rank)
                .map(|node| (node.name.split(' ').count(), node.tax_id)).collect();
            same_rank.sort_by_key(|(count, _)| std::cmp::Reverse(*count));
            debug!("Found tax_id={} {}", matched_nodes[0].tax_id, matched_nodes[0].name);
            return Some(same_rank[0].1);
        }
        None
    }

    /// Access the [`Taxonomy`] object used by this [`TaxonomyMatcher`]
    pub fn taxonomy(&self) -> &Taxonomy { &self.taxonomy }

    fn replace_punctuation_with_spaces(text: &str) -> String {
        let ret = text.chars()
            .map(|c| if c.is_ascii_punctuation() { ' ' } else { c })
            .collect::<String>();
        format!(" {} ", ret)
    }
}