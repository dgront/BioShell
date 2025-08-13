use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use pyo3::types::{PyAny, PyAnyMethods};

use std::sync::Arc;
use once_cell::sync::OnceCell;

use bioshell_taxonomy::{Node, Taxonomy, Rank, TaxonomyMatcher};

use crate::taxonomy::PyRank;

#[pyclass(name = "Node")]
pub struct PyNode(pub(crate) Node);

#[pymethods]
impl PyNode {
    /// Unique taxonomy identifier.
    #[getter]
    fn tax_id(&self) -> u32 { self.0.tax_id }

    /// Parent tax_id in the taxonomy.
    #[getter]
    fn parent_tax_id(&self) -> u32 { self.0.parent_tax_id }

    /// Scientific name.
    #[getter]
    fn name(&self) -> &str { &self.0.name }

    /// Taxonomic rank.
    #[getter]
    fn rank(&self) -> PyRank { PyRank { inner: self.0.rank } }

    fn __repr__(&self) -> String {
        format!(
            "Node(tax_id={}, parent_tax_id={}, rank={}, name='{}')",
            self.0.tax_id,
            self.0.parent_tax_id,
            self.0.rank.to_string(),
            self.0.name
        )
    }
}

impl From<&Node> for PyNode {
    fn from(n: &Node) -> Self {
        let nn = n.clone();
        PyNode(n.clone())
    }
}

/// Python wrapper for the `Taxonomy` struct.
#[pyclass(name = "Taxonomy")]
pub struct PyTaxonomy {
    pub(crate) taxonomy: Arc<Taxonomy>,
    matcher: OnceCell<TaxonomyMatcher>,
}

#[pymethods]
impl PyTaxonomy {
    /// Load taxonomy from a `.tar.gz` archive.
    #[staticmethod]
    pub fn load_from_tar_gz(path: &str) -> PyResult<Self> {
        Taxonomy::load_from_tar_gz(path)
            .map(Arc::new)
            .map(|taxonomy| Self {
                taxonomy,
                matcher: OnceCell::new(),
            })
            .map_err(|e| PyValueError::new_err(format!("Failed to load taxonomy: {}", e)))
    }

    /// Download taxonomy data from NCBI into the specified path.
    #[staticmethod]
    pub fn download_from_ncbi(output_path: &str) -> PyResult<()> {
        Taxonomy::download_from_ncbi(output_path)
            .map_err(|e| PyValueError::new_err(format!("Download failed: {}", e)))
    }

    /// Return a node by its taxid.
    pub fn node(&self, taxid: u32) -> Option<PyNode> {
        self.taxonomy.node(taxid).map(PyNode::from)
    }

    /// Return all nodes as a list.
    pub fn nodes(&self) -> Vec<PyNode> {
        self.taxonomy.nodes().map(PyNode::from).collect()
    }

    /// Return all names associated with a taxid.
    pub fn names(&self, taxid: u32) -> Vec<String> {
        self.taxonomy.names(taxid).cloned().collect()
    }

    /// Return the full lineage (list of nodes) for a taxid.
    pub fn lineage(&self, taxid: u32) -> Vec<PyNode> {
        self.taxonomy.lineage(taxid).into_iter().map(PyNode::from).collect()
    }

    /// Look up a taxid by scientific name.
    pub fn taxid(&self, name: &str) -> Option<u32> {
        self.taxonomy.taxid(name)
    }

    /// Return the node in the lineage with the given rank.
    pub fn rank<'py>(&self, taxid: u32, rank: &Bound<'py, PyAny>) -> PyResult<Option<PyNode>> {
        if let Ok(py_rank) = rank.extract::<PyRef<PyRank>>() {
            Ok(self.taxonomy.rank(taxid, py_rank.inner).map(PyNode::from))
        } else if let Ok(s) = rank.extract::<&str>() {
            Ok(self.taxonomy.rank(taxid, Rank::from_str(s)).map(PyNode::from))
        } else {
            Err(PyValueError::new_err("Expected a Rank or str as the second argument"))
        }
    }

    pub fn classification(&self, tax_id: u32) -> PyResult<Vec<Option<PyNode>>> {

        let mut output = vec![];
        for rnk in self.taxonomy.classification(tax_id).into_iter() {
            if let Some(n) = rnk {
                output.push(Some(PyNode::from(n)))  ;
            } else {
                output.push(None);
            }
        }

        Ok(output)
    }

    /// Find a species information in a given string.
    ///
    pub fn find(&self, description: &str) -> Option<u32> {
        if let Some(matcher_ref) = self.matcher.get() {
            matcher_ref.find(description)
        } else {
            let matcher = TaxonomyMatcher::from_arc(Arc::clone(&self.taxonomy))
                .unwrap_or_else(|e| panic!("Failed to initialize matcher: {}", e));

            self.matcher
                .set(matcher)
                .unwrap_or_else(|_| panic!("Matcher was already initialized unexpectedly"));

            self.matcher.get().unwrap().find(description)
        }
    }
}

