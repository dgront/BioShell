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
    /// be more verbose and log program actions on the screen
    #[clap(short, long, short='v')]
    verbose: bool
}

fn main() {
    let args = Args::parse();
    unsafe {
        if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
        if args.verbose { env::set_var("RUST_LOG", "debug"); }
    }
    env_logger::init();

}