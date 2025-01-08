use nalgebra::*;

pub fn run_hpp_attack() {
    // runs the HPP attack against Hawk
    //
    // Inputs: t digital signatures each signed with a secret key B
    // Outputs: Columns of B, or failure

    // In real life attack one would obtain the signatures as 
    // sig = enc(s1) where s1 = (h1 - w1) / 2
    // and one would need to recompute s0 on the signer/attacker side.
    // However for this implementation we simply generate entire signatures on signer end 
    // and directly return the entire signature w which we use to deduce information about B

    // STEP 0: Generate samples
    // We need to generate a lot of samples to run the attack on

    // STEP 1: estimate covariance matrix. This step requires a lot of samples,
    // so hopefully we can employ some sort of trick like the HPP against NTRU to reduce
    // number of signatures needed
    
    // STEP 2: conversion from hidden parallelepiped to hidden hypercube.  
    // in this step we need a covariance matrix estimation from step 1. The better
    // the estimation in step two, the better conversion estimation we can do here.
    //
    // Given matrix G, we want to compute L s.t. L^t L = G^-1, and thereafter 
    // multiply our signature samples on the right with this L 
    // We use the Nalgebra crate for representations of matrices and for procedures such
    // as Cholesky decomposition
    
    // STEP 3: Gradient Descent:
    // The final step is to do gradient descent on our (converted) samples to minimize the 
    // fourth moment, and consequently reveal a row/column from +/- B
}

fn generate_samples(t: usize, degree: usize) -> Vec<i16> {
    vec![]
}

fn estimate_covariance_matrix(samples: Vec<i16>) -> DMatrix<i16> {
    let g: DMatrix<i16> = DMatrix::identity(100, 100);
    g
}

fn hypercube_transformation(samples: Vec<i16>, g: DMatrix<i16>) -> Vec<i16> {
    samples
}

fn gradient_descent(samples: Vec<i16>) {

}
