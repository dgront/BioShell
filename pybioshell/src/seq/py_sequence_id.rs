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

    /// Create a new `TrEMBL` id from a given string.
    ///
    /// `TrEMBL` id actually is the `UniProtKB` id with the `tr` prefix to indicate that it is a TrEMBL entry.
    #[staticmethod]
    pub fn trembl(value: String) -> Self {
        Self { inner: SeqId::TrEmbl(value) }
    }

    /// Create a new `UniProtKB` id from a given string.
    #[staticmethod]
    pub fn uniprot_id(value: String) -> Self {
        Self { inner: SeqId::UniProtKB(value) }
    }

    /// Create a new `UniRef` id from a given string.
    ///
    /// UniRef id is a cluster of UniProtKB entries, and the id is usually in the form of `UniRef100_P12345`,
    /// where `P12345` is a UniProtKB id and UniRef100 shows the clustering was based on 100% sequence identity cutoff.
    #[staticmethod]
    pub fn uniref(value: String) -> Self {
        Self { inner: SeqId::UniRef(value) }
    }

    /// Create a new `UniParc` id from a given string.
    ///
    /// UniParc is the unique identifier assigned to a distinct
    /// protein sequence. It consists of the characters “UPI” followed
    /// by 10 hexadecimal characters (0–9, A–F). A UPI is stable across releases.
    #[staticmethod]
    pub fn uniparc(value: String) -> Self {
        Self { inner: SeqId::UniParc(value) }
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

    /// Create a new `KEGG` id from a given string.
    #[staticmethod]
    pub fn kegg(value: String) -> Self {
        Self { inner: SeqId::KEGG(value) }
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

    /// Create a new `CypId` id from a given string.
    #[staticmethod]
    pub fn cypid(value: String) -> Self {
        Self { inner: SeqId::CypId(value) }
    }

    /// Create a new `DDBJ` id from a given `value` string.
    #[staticmethod]
    pub fn ddbj(value: String) -> Self {
        Self { inner: SeqId::DDBJ(value) }
    }


    /// Return the kind of this sequence id (e.g., "PDB", "SwissProt", etc.).
    ///
    /// Currently, bioshell recognizes the following kinds of sequence identifiers:
    /// "PDB", "SwissProt", "UniProtKB", "UniRef", "RefSeq", "GenBank", "Ensembl", "TrEMBL", "NCBIGI",
    /// "UniProtEntry", "DDBJ", "TaxId", "Organism", and "CypId". There is also a "Default" kind
    /// for identifiers that do not match any of the above.
    #[getter]
    pub fn kind(&self) -> &'static str {
        match &self.inner {
            SeqId::PDB(_) => "PDB",
            SeqId::SwissProt(_) => "SwissProt",
            SeqId::UniProtKB(_) => "UniProtKB",
            SeqId::UniParc(_) => "UniParc",
            SeqId::UniRef(_) => "UniRef",
            SeqId::RefSeq(_) => "RefSeq",
            SeqId::GenBank(_) => "GenBank",
            SeqId::Ensembl(_) => "Ensembl",
            SeqId::TrEmbl(_) => "TrEMBL",
            SeqId::NCBIGI(_) => "NCBIGI",
            SeqId::KEGG(_) => "KEGG",
            SeqId::TaxId(_) => "TaxId",
            SeqId::Organism(_) => "Organism",
            SeqId::UniProtEntry(_) => "UniProtEntry",
            SeqId::DDBJ(_) => "DDBJ",
            SeqId::CypId(_) => "CypId",
            SeqId::Default(_) => "Default",
        }
    }

    /// Return the actual identifier string, e.g., "2gb1" for a PDB id, "P12345" for a SwissProt id, etc.
    #[getter]
    pub fn value(&self) -> &str {
        match &self.inner {
            SeqId::PDB(v)
            | SeqId::SwissProt(v)
            | SeqId::UniProtKB(v)
            | SeqId::UniParc(v)
            | SeqId::UniRef(v)
            | SeqId::RefSeq(v)
            | SeqId::GenBank(v)
            | SeqId::Ensembl(v)
            | SeqId::TrEmbl(v)
            | SeqId::NCBIGI(v)
            | SeqId::TaxId(v)
            | SeqId::KEGG(v)
            | SeqId::Organism(v)
            | SeqId::UniProtEntry(v)
            | SeqId::DDBJ(v)
            | SeqId::CypId(v)
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

#[doc = include_str!("../../docs/seq/parse_sequence_id.rst")]
#[pyfunction]
pub fn parse_sequence_id<'py>(py: Python<'py>,  description: &str) -> PyResult<Bound<'py, PyList>> {

    let ids: SeqIdList = bioshell_seq::sequence::parse_sequence_id(description);

    let py_ids: Vec<PySeqId> = ids.into_iter().map(|inner| PySeqId { inner }).collect();

    PyList::new(py, py_ids)
}