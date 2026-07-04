use pyo3::prelude::*;
use pyo3::types::PyList;

use bioshell_seq::sequence::{SeqId, SeqIdList};


#[pyclass(name = "SeqId")]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PySeqId {
    inner: SeqId,
}

#[pymethods]
impl PySeqId {

    #[doc = include_str!("../../docs/seq/seq_id.pdb.rst")]
    #[staticmethod]
    pub fn pdb(value: String) -> Self {
        Self { inner: SeqId::PDB(value) }
    }

    #[doc = include_str!("../../docs/seq/seq_id.swissprot.rst")]
    #[staticmethod]
    pub fn swissprot(value: String) -> Self {
        Self { inner: SeqId::SwissProt(value) }
    }

    /// Create a new `UniProt` id from a given string.
    #[staticmethod]
    pub fn uniprot_id(value: String) -> Self {
        Self { inner: SeqId::UniProtID(value) }
    }

    /// Create a new `UniRef` id from a given string.
    #[staticmethod]
    pub fn uniref(value: String) -> Self {
        Self { inner: SeqId::UniRef(value) }
    }

    /// Create a new `RefSeq` id from a given string.
    #[staticmethod]
    pub fn refseq(value: String) -> Self {
        Self { inner: SeqId::RefSeq(value) }
    }

    /// Create a new `GeneBank` id from a given string.
    #[staticmethod]
    pub fn genbank(value: String) -> Self {
        Self { inner: SeqId::GenBank(value) }
    }

    /// Create a new `Ensembl` id from a given string.
    #[staticmethod]
    pub fn ensembl(value: String) -> Self {
        Self { inner: SeqId::Ensembl(value) }
    }

    /// Create a new `TrEMBL` id from a given string.
    #[staticmethod]
    pub fn trembl(value: String) -> Self {
        Self { inner: SeqId::TrEmbl(value) }
    }

    /// Create a new `NCBIgi` id from a given string.
    #[staticmethod]
    pub fn ncbi_gi(value: String) -> Self {
        Self { inner: SeqId::NCBIGI(value) }
    }

    /// Create a new `TaxID` id from a given string.
    #[staticmethod]
    pub fn taxid(value: String) -> Self {
        Self { inner: SeqId::TaxId(value) }
    }

    /// Create a new id without specifying which database it refers to
    #[staticmethod]
    pub fn default(value: String) -> Self {
        Self { inner: SeqId::Default(value) }
    }

    /// Return the kind of this sequence id (e.g., "PDB", "SwissProt", etc.)
    #[getter]
    pub fn kind(&self) -> &'static str {
        match &self.inner {
            SeqId::PDB(_) => "PDB",
            SeqId::SwissProt(_) => "SwissProt",
            SeqId::UniProtID(_) => "UniProtID",
            SeqId::UniRef(_) => "UniRef",
            SeqId::RefSeq(_) => "RefSeq",
            SeqId::GenBank(_) => "GenBank",
            SeqId::Ensembl(_) => "Ensembl",
            SeqId::TrEmbl(_) => "TrEmbl",
            SeqId::NCBIGI(_) => "NCBIGI",
            SeqId::TaxId(_) => "TaxId",
            SeqId::Default(_) => "Default",
        }
    }

    /// Return the actual identifier string, e.g., "2gb1" for a PDB id, "P12345" for a SwissProt id, etc.
    #[getter]
    pub fn value(&self) -> &str {
        match &self.inner {
            SeqId::PDB(v)
            | SeqId::SwissProt(v)
            | SeqId::UniProtID(v)
            | SeqId::UniRef(v)
            | SeqId::RefSeq(v)
            | SeqId::GenBank(v)
            | SeqId::Ensembl(v)
            | SeqId::TrEmbl(v)
            | SeqId::NCBIGI(v)
            | SeqId::TaxId(v)
            | SeqId::Default(v) => v,
        }
    }

    /// Return a tuple of (kind, value) for this sequence id.
    pub fn as_tuple(&self) -> (&'static str, &str) {
        (self.kind(), self.value())
    }

    fn __repr__(&self) -> String { format!("{}", &self.inner) }

    fn __str__(&self) -> String { format!("{}", &self.inner) }
}

#[pyfunction]
pub fn parse_sequence_id<'py>(py: Python<'py>,  description: &str) -> PyResult<Bound<'py, PyList>> {

    let ids: SeqIdList = bioshell_seq::sequence::parse_sequence_id(description);

    let py_ids: Vec<PySeqId> = ids.into_iter().map(|inner| PySeqId { inner }).collect();

    PyList::new(py, py_ids)
}