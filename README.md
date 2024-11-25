Master project for Eirik Skjerve, Secure and Reliable communication, Department of Informatics, University of Bergen
---
### Brief description
This project will revolve around implementation and analysis of [HAWK](https://hawk-sign.info/), a digital signature scheme submitted in the NIST-competition for a new digital signature standard.
Hawk is based upon the Lattice Isomorphism Problem, a problem believed to be hard to solve even for quantum computers.

### How to run Hawk benchmarks
1. Make sure Rust and Cargo is installed
2. Clone this repo
3. cd hawk_rust
4. cargo run (optional parameter "--release" for optimized compilation)
<br/>
Average timing result of key generation and 100 signature generations and verifications for 
each degree 256, 512, and 1024

### Roadmap
- I have implemented Hawk in the [Rust](https://www.rust-lang.org/) programming language.
- I wish later to implement a cryptanalytic attack based on [Learning a Parallelepiped](https://cims.nyu.edu/~regev/papers/gghattack.pdf) to attempt to recover (parts of) a secret key given enough signatures.
