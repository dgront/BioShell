use std::str::FromStr;
use bioshell_pdb::{Deposit, list_ligands_in_deposit};

fn deposit_title(deposit: &Deposit) -> String { deposit.title.clone().unwrap_or("".to_string()) }
fn deposit_resolution(deposit: &Deposit) -> String {
    if let Some(res) = deposit.resolution {res.to_string()}
    else { "".to_string() }
}

fn deposit_methods(deposit: &Deposit) -> String {
    deposit.methods.iter().map(|m|m.to_string()).collect::<Vec<_>>().join(",")
}

fn deposit_ligands(deposit: &Deposit) -> String {
    list_ligands_in_deposit(deposit).iter()
        .map(|l| l.code3.clone())
        .collect::<Vec<String>>().join(" ")
}
fn deposit_classification(deposit: &Deposit) -> String { deposit.classification.clone().unwrap_or("".to_string()) }
fn deposit_id(deposit: &Deposit) -> String { deposit.id_code.clone() }
fn deposit_r_factor(deposit: &Deposit) -> String {
    if let Some(r_fact) = deposit.r_factor { r_fact.to_string() }
     else { "".to_string() }
}
fn deposit_keywords(deposit: &Deposit) -> String { deposit.keywords.join(",") }
fn space_group(deposit: &Deposit) -> String { deposit.unit_cell.as_ref().unwrap().space_group.clone() }

const DEPOSIT_INFO_FEATURES: [&str; 10] = ["id", "keywords", "title", "resolution", "methods",
                                "ligands", "rfactor", "spacegroup", "unitcell", "classification"];

pub(crate) fn get_deposit_info<'a>(dep: &'a Deposit, deposit_info_names: &'a Vec<String>) -> Vec<(&'a str, String)> {

    let tokens = if deposit_info_names.is_empty() {
        &DEPOSIT_INFO_FEATURES.into_iter()
            .map(String::from)
            .collect()
    } else { deposit_info_names };

    let mut out = vec![];
    for name in tokens {

        match DepositInfo::from_str(name) {
            Ok(info) => {
                // Match on each DepositInfo variant to access and call the function
                let result = match info {
                    DepositInfo::Id((token_name, func)) => (token_name, func(dep)),
                    DepositInfo::Keywords((token_name, func)) => (token_name, func(dep)),
                    DepositInfo::Title((token_name, func)) => (token_name, func(dep)),
                    DepositInfo::Resolution((token_name, func)) => (token_name, func(dep)),
                    DepositInfo::Methods((token_name, func)) => (token_name, func(dep)),
                    DepositInfo::Ligands((token_name, func)) => (token_name, func(dep)),
                    DepositInfo::RFactor((token_name, func)) => (token_name, func(dep)),
                    DepositInfo::Classification((token_name, func)) => (token_name, func(dep)),
                    DepositInfo::SpaceGroup((token_name, func)) => (token_name, func(dep)),
                    DepositInfo::UnitCell((token_name, func)) => (token_name, func(dep)),
                };
                out.push(result);
            }
            Err(_) => { eprintln!("Error: Unknown deposit property: {}; allowed keywords are: {:?}", name, DEPOSIT_INFO_FEATURES); }
        }
    }

    return out;
}

#[derive(Debug, Clone, Copy)]
pub(crate) enum DepositInfo {
    Id((&'static str, fn(&Deposit) -> String)),
    Keywords((&'static str, fn(&Deposit) -> String)),
    Title((&'static str, fn(&Deposit) -> String)),
    Resolution((&'static str, fn(&Deposit) -> String)),
    Methods((&'static str, fn(&Deposit) -> String)),
    Ligands((&'static str, fn(&Deposit) -> String)),
    RFactor((&'static str, fn(&Deposit) -> String)),
    Classification((&'static str, fn(&Deposit) -> String)),
    SpaceGroup((&'static str, fn(&Deposit) -> String)),
    UnitCell((&'static str, fn(&Deposit) -> String)),
}

impl FromStr for DepositInfo {
    type Err = &'static str;

    fn from_str(input: &str) -> Result<Self, Self::Err> {
        match input.to_lowercase().as_str() {
            "id" => Ok(DepositInfo::Id(("id", deposit_id))),
            "keywords" => Ok(DepositInfo::Keywords(("keywords", deposit_keywords))),
            "title" => Ok(DepositInfo::Title(("title", deposit_title))),
            "resolution" => Ok(DepositInfo::Resolution(("resolution", deposit_resolution))),
            "methods" => Ok(DepositInfo::Methods(("methods", deposit_methods))),
            "ligands" => Ok(DepositInfo::Ligands(("ligands", deposit_ligands))),
            "rfactor" => Ok(DepositInfo::RFactor(("rfactor", deposit_r_factor))),
            "spacegroup" => Ok(DepositInfo::SpaceGroup(("space group", deposit_ligands))),
            "unitcell" => Ok(DepositInfo::UnitCell(("unit cell", deposit_r_factor))),
            "classification" => Ok(DepositInfo::Classification(("classification", deposit_classification))),
            _ => Err("Error: Unknown deposit property"),
        }
    }
}
