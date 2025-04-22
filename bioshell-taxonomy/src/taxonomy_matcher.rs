use std::error::Error;
use std::time::Instant;
use aho_corasick::AhoCorasick;
use log::{debug, info};
use crate::{Node, Taxonomy};

pub struct TaxonomyMatcher<'a> {
    taxonomy: &'a Taxonomy,
    matcher: AhoCorasick,
    patterns: Vec<String>, // store names for mapping back
}


impl<'a> TaxonomyMatcher<'a> {
    pub fn new(taxonomy: &'a Taxonomy) -> Result<Self, Box<dyn Error>> {

        let start = Instant::now();
        let patterns = taxonomy.nodes().map(|n| n.name.clone()).collect::<Vec<String>>();
        let matcher = AhoCorasick::new(&patterns)?;
        info!("Taxonomy matcher constructed in {:.2?}", start.elapsed());
        Ok(TaxonomyMatcher { taxonomy, matcher, patterns })
    }

    pub fn find(&self, description: &str) -> Option<u32> {

        let desc = description.replace("_", " ");
        let mut matches = self.matcher.find_overlapping_iter(&desc);
        let mut matched_nodes: Vec<&Node> = vec![];
        while let Some(m) = matches.next() {
            let detected_name = &self.patterns[m.pattern()];
            let taxid = self.taxonomy.taxid(detected_name).unwrap();
            matched_nodes.push(self.taxonomy.node(taxid).unwrap());
        }
        if matched_nodes.len() > 0 {
            matched_nodes.sort_by_key(|node| node.rank);
            debug!("Found tax_id={} {}", matched_nodes[0].tax_id, matched_nodes[0].name);
            return Some(matched_nodes[0].tax_id);
        }
        None
    }
}