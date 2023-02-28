# BioShell v.4


## What is BioShell?
BioShell project, started in 2006) has been a command line toolkit for  structural bioinformatics.
Over time it has also provided Coarse Grained simulations of proteins

## WARNING
The fourth version of the BioShell package is still in the _very_ early stage of development.
Many important features are still missing and APIs can and will change. 
Documentation is sparse. If you are looking for a stable package, ypu might consider using
the BioShell v.3 package. 

## Project structure
Curretnly the BioShell v.4 project has been divided into the following crates:

 - **bioshell-core**: core BioShell's functionality
 - **bioshell-sim**: base traits (interfaces) defining molecular simulation engine 
 - **bioshell-montecarlo**: Monte Carlo sampling schemes
 - **bioshell-statistics**: statistical utilities used in bioinformatics applications as well as during simulations
 - **bioshell-numerical**: data structures (such as a 3D vector) and numerical utilities that operate on them 
 - **bioshell-clustering**: clustering methods, such as K-means and OPTICS
 - **bioshell-cartesians**: provides structs necessary for modeling Cartesian systems, such as polymers, argon fluid, etc.
 - **bioshell-ff**: provides molecular force fields

## Building
You need to install `rust` toolchain to compile the package.
- to compile all libaries and executables:
```bash
cargo build --release
```
 - to compile an example:
```bash
cargo build --example polymer_perm --release
```
- to compile all examples:
```bash
cargo build --examples --release
```

 - to build documentation:
```
RUSTDOCFLAGS="--html-in-header katex-header.html" cargo doc --no-deps --open
```