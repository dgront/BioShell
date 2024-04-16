use rand::rngs::SmallRng;
use rand::SeedableRng;
use surpass::{BoxDefinition, SurpassAlphaSystem};
#[test]
fn surpass_alpha_from_fasta() {

    let mut rnd = SmallRng::from_entropy();
    let sec_str = vec!["CEEEEEECCCCCCEEEEEECCHHHHHHHHHHHHHHHCCCCCEEEEECCCCEEEEEC".to_string()];
    let system = SurpassAlphaSystem::by_secondary_structure(&sec_str, BoxDefinition::Density {density:0.001}, &mut rnd);
}