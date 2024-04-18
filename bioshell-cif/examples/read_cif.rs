use std::env;
use bioshell_cif::{read_cif_file};

fn main() {
    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    env_logger::init();

    let args: Vec<String> = env::args().collect();
    let cif_content = read_cif_file(&args[1]).unwrap();
    println!("Number of data blocks: {}", cif_content.len());

    for block in &cif_content {
        println!("data block: {}: {:4} entries and {:4}  loops",
                 block.name(), block.data_items().len(), block.loop_blocks().count());
    }
}
