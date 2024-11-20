Master project for Eirik Skjerve, Secure and Reliable communication, Department of Informatics, University of Bergen
---
### Brief description
This project will revolve around implementation and analysis of a selection of the digital signature schemes submitted in the NIST-competition for a new digital signature standard.
We are currently interested in [HAWK](https://hawk-sign.info/), which is based on the Lattice Isomorphism Problem, a problem believed to be hard to solve even for quantum computers. 

### How to run Hawk benchmarks
1. Make sure Rust and Cargo is installed
2. Clone this repo
3. cd hawk_rust
4. cargo run (optional parameter "--release" for optimized compilation)
Average timing result of key generation and 100 signature generations and verifications for 
each degree 256, 512, and 1024

### Roadmap
- I am currently working on implementing the HAWK scheme in the [Rust](https://www.rust-lang.org/) programming language.
- I wish later to implement a cryptanalytic attack based on [HPP](https://cims.nyu.edu/~regev/papers/gghattack.pdf) (Hidden Parallelepiped Problem) to attempt to recover (parts of) a secret key given enough signatures generated by the secret key.
