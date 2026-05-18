use std::fmt;
use std::str::FromStr;
use crate::ChemErrors;

macro_rules! define_elements {
    (
        $(
            $variant:ident = $z:literal, $symbol:literal, $name:literal, $mass:literal, $valence:expr;
        )+
    ) => {
        /// Defines atomic elements and their basic properties
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
        #[repr(u8)]
        pub enum Element {
            $(
                $variant = $z,
            )+
        }

        impl Element {
            pub const ALL: &'static [Element] = &[
                $(
                    Element::$variant,
                )+
            ];

            /// Atomic number of this element
            pub fn atomic_number(self) -> u8 {
                self as u8
            }

            /// Chemical symbol of this element
            pub fn symbol(self) -> &'static str {
                match self {
                    $(
                        Element::$variant => $symbol,
                    )+
                }
            }

            /// THe name of this element in plain English
            pub fn name(self) -> &'static str {
                match self {
                    $(
                        Element::$variant => $name,
                    )+
                }
            }

            /// Atomic (molar) mass of this element
            pub fn mass(self) -> f64 {
                match self {
                    $(
                        Element::$variant => $mass,
                    )+
                }
            }

            pub fn default_valence(self) -> Option<u8> {
                match self {
                    $(
                        Element::$variant => $valence,
                    )+
                }
            }

            pub fn from_atomic_number(z: u8) -> Result<Self, ChemErrors> {
                match z {
                    $(
                        $z => Ok(Element::$variant),
                    )+
                    _ => Err(ChemErrors::InvalidAtomicNumber(z)),
                }
            }
        }

        impl FromStr for Element {
            type Err = ChemErrors;

            fn from_str(s: &str) -> Result<Self, Self::Err> {
                match s.trim() {
                    $(
                        $symbol => Ok(Element::$variant),
                    )+
                    other => Err(ChemErrors::UnknownElement(other.to_string())),
                }
            }
        }

        impl fmt::Display for Element {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_str(self.symbol())
            }
        }
    };
}

define_elements! {
    H  = 1,  "H",  "hydrogen",     1.008,    Some(1);
    He = 2,  "He", "helium",       4.002602, None;
    Li = 3,  "Li", "lithium",      6.94,     Some(1);
    Be = 4,  "Be", "beryllium",    9.012183, Some(2);
    B  = 5,  "B",  "boron",        10.81,    Some(3);
    C  = 6,  "C",  "carbon",       12.011,   Some(4);
    N  = 7,  "N",  "nitrogen",     14.007,   Some(3);
    O  = 8,  "O",  "oxygen",       15.999,   Some(2);
    F  = 9,  "F",  "fluorine",     18.998403163, Some(1);
    Ne = 10, "Ne", "neon",         20.1797,  None;

    Na = 11, "Na", "sodium",       22.98976928, Some(1);
    Mg = 12, "Mg", "magnesium",    24.305,      Some(2);
    Al = 13, "Al", "aluminium",    26.9815385,  Some(3);
    Si = 14, "Si", "silicon",      28.085,      Some(4);
    P  = 15, "P",  "phosphorus",   30.973761998, Some(3);
    S  = 16, "S",  "sulfur",       32.06,       Some(2);
    Cl = 17, "Cl", "chlorine",     35.45,       Some(1);
    Ar = 18, "Ar", "argon",        39.948,      None;

    K  = 19, "K",  "potassium",    39.0983, Some(1);
    Ca = 20, "Ca", "calcium",      40.078,  Some(2);
    Sc = 21, "Sc", "scandium",     44.955908, None;
    Ti = 22, "Ti", "titanium",     47.867,    None;
    V  = 23, "V",  "vanadium",     50.9415,   None;
    Cr = 24, "Cr", "chromium",     51.9961,   None;
    Mn = 25, "Mn", "manganese",    54.938044, None;
    Fe = 26, "Fe", "iron",         55.845,    None;
    Co = 27, "Co", "cobalt",       58.933194, None;
    Ni = 28, "Ni", "nickel",       58.6934,   None;
    Cu = 29, "Cu", "copper",       63.546,    None;
    Zn = 30, "Zn", "zinc",         65.38,     Some(2);

    Ga = 31, "Ga", "gallium",      69.723,  Some(3);
    Ge = 32, "Ge", "germanium",    72.630,  Some(4);
    As = 33, "As", "arsenic",      74.921595, Some(3);
    Se = 34, "Se", "selenium",     78.971,  Some(2);
    Br = 35, "Br", "bromine",      79.904,  Some(1);
    Kr = 36, "Kr", "krypton",      83.798,  None;

    Rb = 37, "Rb", "rubidium",     85.4678, Some(1);
    Sr = 38, "Sr", "strontium",    87.62,   Some(2);
    Y  = 39, "Y",  "yttrium",      88.90584, None;
    Zr = 40, "Zr", "zirconium",    91.224,   None;
    Nb = 41, "Nb", "niobium",      92.90637, None;
    Mo = 42, "Mo", "molybdenum",   95.95,    None;
    Tc = 43, "Tc", "technetium",   98.0,     None;
    Ru = 44, "Ru", "ruthenium",    101.07,   None;
    Rh = 45, "Rh", "rhodium",      102.90550, None;
    Pd = 46, "Pd", "palladium",    106.42,   None;
    Ag = 47, "Ag", "silver",       107.8682, Some(1);
    Cd = 48, "Cd", "cadmium",      112.414,  Some(2);

    In = 49, "In", "indium",       114.818, Some(3);
    Sn = 50, "Sn", "tin",          118.710, Some(4);
    Sb = 51, "Sb", "antimony",     121.760, Some(3);
    Te = 52, "Te", "tellurium",    127.60,  Some(2);
    I  = 53, "I",  "iodine",       126.90447, Some(1);
    Xe = 54, "Xe", "xenon",        131.293, None;
}