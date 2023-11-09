# BioShell v.4


## What is BioShell?
BioShell project, started in 2006, has been a command line toolkit for  structural bioinformatics.
Over time, it has also provided Coarse Grained simulations of proteins

## WARNING
The fourth version of the BioShell package is still in the _very_ early stage of development.
Many important features are still missing and APIs can and will change. 
Documentation is sparse. If you are looking for a stable package, ypu might consider using
the BioShell v.3 package. 

## Project structure
Curretnly the BioShell v.4 project has been divided into the following crates:

 - **bioshell-seq**: works with biological sequences
 - **bioshell-pdb**: a library to work on with crystallographic Protein DataBank files. It can parse PDB files and perform calculations on macromolecular structures.
 - **bioshell-cif**: reads CIF files
 - **bioshell-statistics**: statistical utilities 
 - **bioshell-clustering**: clustering methods, such as K-means and OPTICS
 - **bioshell-io**: I/O utilities

## Building
You need to install `rust` toolchain to compile the package. You can:
 
 - compile the whole project, i.e. all its libraries (called crates) by executing the following command in the root folder of the project: 
```bash
cargo build --release
```
This will compile also all the executables (currently the ``surpass_alfa`` application)
- compile all examples:
```bash
cargo build --examples --release
```
- compile a selected example:
```bash
cargo build --example polymer_perm --release
```
All the above commands compile the *release* (i.e. optimized) version. To obtain the *debug* build, just remove the ``--release`` flag.

Documentation has to be built separately for each crate. Enter the crate folder (e.g. ``bioshell-statistics``) and execute the following:
```
RUSTDOCFLAGS="--html-in-header katex-header.html" cargo doc --no-deps --open
```

To run all tests (for all the crates of the project), run the:
```
cargo test
```
command in the root folder. Alternatively, you can test a single crate by running that command in a respective subfolder.