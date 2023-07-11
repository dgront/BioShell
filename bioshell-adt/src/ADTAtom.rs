use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::convert::TryFrom;

#[derive(Debug, PartialEq, Eq, Hash)]
pub struct ADTAtom {
    pub symbol: AtomicSymbolEnum,
    pub position: Option<AtomicPositionEnum>,
    pub branch: Option<i32>,
    pub protein_name: String,
    pub chain_name: String,
    pub residue_no: i32,
    pub atom_no: i32,
    pub name: AtomicNameEnum,
    pub alt_loc: String,
    pub radius: f64,
    pub atomic_number: f64,
    pub atomic_mass: f64,
    pub color_rgb: Color,
    pub coordinate: Point3d,
    pub valency: Vec<i32>,
}

impl ADTAtom {
    pub fn new() -> Self {
        Self {
            symbol: AtomicSymbolEnum::default(),
            position: None,
            branch: None,
            protein_name: String::new(),
            chain_name: String::new(),
            residue_no: 0,
            atom_no: 0,
            name: AtomicNameEnum::default(),
            alt_loc: String::new(),
            radius: 0.0,
            atomic_number: 0.0,
            atomic_mass: 0.0,
            color_rgb: Color::default(),
            coordinate: Point3d::default(),
            valency: Vec::new(),
        }
    }

    pub fn from_text(&mut self, atom_line: &str) {
        let elements = PdbLineParser::parse_atom(atom_line);

        self.atom_no = i32::from_str_radix(elements[1].trim(), 10).unwrap();
        self.name = AtomicNameEnum::try_from(elements[2].trim()).unwrap();

        self.chain_name = elements[5].trim().to_owned();
        self.residue_no = i32::from_str_radix(elements[6].trim(), 10).unwrap();

        self.coordinate = Point3d::new(
            f64::from_str(elements[8].trim()).unwrap(),
            f64::from_str(elements[9].trim()).unwrap(),
            f64::from_str(elements[10].trim()).unwrap(),
        );
    }

    pub fn new_with_params(name: &str, radius: f64, rgb: Color) -> Self {
        Self {
            name: AtomicNameEnum::try_from(name).unwrap(),
            radius,
            color_rgb: rgb,
            ..Default::default()
        }
    }

    pub fn set_dic_element(&mut self, at_element: &AtElement) {
        self.name = at_element.name;
        self.symbol = at_element.symbol;
        self.valency = at_element.valency.clone();
        self.radius = at_element.radius;
        self.atomic_number = at_element.atomic_number;
        self.atomic_mass = at_element.atomic_mass;
        self.color_rgb = at_element.color_rgb;
    }

    pub fn set_pdb(&mut self, pdb_atom: &PdbAtom) {
        self.protein_name = pdb_atom.protein_name.clone();
        self.residue_no = pdb_atom.residue_no.unwrap();
        self.alt_loc = pdb_atom.alt_loc_indicator.clone();
        self.atom_no = pdb_atom.atom_serial_no.unwrap();

        self.symbol = AtomicSymbolEnum::try_from(pdb_atom.atom_symbol.as_str()).unwrap();
        self.position = match pdb_atom.atom_position.as_ref() {
            Some(pos) => Some(AtomicPositionEnum::try_from(pos.as_str()).unwrap()),
            None => None,
        };
        self.branch = pdb_atom.atom_no_in_the_branch;

        self.chain_name = pdb_atom.chain_name.clone();
        self.residue_no = pdb_atom.residue_no.unwrap();
        self.coordinate = pdb_atom.coordinate;
    }
}

impl Hash for ADTAtom {
    fn hash<H: Hasher>(&self, state: &mut H) {
        (
            &self.protein_name,
            &self.chain_name,
            &self.residue_no,
            &self.atom_no,
            &self.name,
            &self.alt_loc,
            &self.radius,
            &self.atomic_number,
            &self.atomic_mass,
        )
            .hash(state);
    }
}

impl PartialEq for ADTAtom {
    fn eq(&self, other: &Self) -> bool {
        self.protein_name == other.protein_name
            && self.chain_name == other.chain_name
            && self.residue_no == other.residue_no
            && self.name == other.name
    }
}

impl Eq for ADTAtom {}

#[derive(Debug, PartialEq, Eq, Hash)]
pub enum AtomicSymbolEnum {
    // Define your enum variants here
}

impl Default for AtomicSymbolEnum {
    fn default() -> Self {
        Self::A // Set a default variant
    }
}