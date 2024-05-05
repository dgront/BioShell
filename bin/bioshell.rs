use clap::Parser;
use std::env;
#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
/// the main BioShell program
/// say bioshell -h to see options
struct Args {
    /// input file name
    #[clap(short, long, short='i', required=true)]
    infile: String,
}

fn main() {
    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    env_logger::init();

    let args = Args::parse();
}