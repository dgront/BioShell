use std::fmt;
use std::str::FromStr;
use crate::ChemErrors;

macro_rules! define_elements {
    (
        $(
            $variant:ident = $z:literal, $symbol:literal, $name:literal;
        )+
    ) => {
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

            pub fn atomic_number(self) -> u8 {
                self as u8
            }

            pub fn symbol(self) -> &'static str {
                match self {
                    $(
                        Element::$variant => $symbol,
                    )+
                }
            }

            pub fn name(self) -> &'static str {
                match self {
                    $(
                        Element::$variant => $name,
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
    H  = 1,  "H",  "hydrogen";
    He = 2,  "He", "helium";
    Li = 3,  "Li", "lithium";
    Be = 4,  "Be", "beryllium";
    B  = 5,  "B",  "boron";
    C  = 6,  "C",  "carbon";
    N  = 7,  "N",  "nitrogen";
    O  = 8,  "O",  "oxygen";
    F  = 9,  "F",  "fluorine";
    Ne = 10, "Ne", "neon";

    Na = 11, "Na", "sodium";
    Mg = 12, "Mg", "magnesium";
    Al = 13, "Al", "aluminium";
    Si = 14, "Si", "silicon";
    P  = 15, "P",  "phosphorus";
    S  = 16, "S",  "sulfur";
    Cl = 17, "Cl", "chlorine";
    Ar = 18, "Ar", "argon";

    K  = 19, "K",  "potassium";
    Ca = 20, "Ca", "calcium";
    Sc = 21, "Sc", "scandium";
    Ti = 22, "Ti", "titanium";
    V  = 23, "V",  "vanadium";
    Cr = 24, "Cr", "chromium";
    Mn = 25, "Mn", "manganese";
    Fe = 26, "Fe", "iron";
    Co = 27, "Co", "cobalt";
    Ni = 28, "Ni", "nickel";
    Cu = 29, "Cu", "copper";
    Zn = 30, "Zn", "zinc";

    Ga = 31, "Ga", "gallium";
    Ge = 32, "Ge", "germanium";
    As = 33, "As", "arsenic";
    Se = 34, "Se", "selenium";
    Br = 35, "Br", "bromine";
    Kr = 36, "Kr", "krypton";
}
