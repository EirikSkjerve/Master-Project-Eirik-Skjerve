use nalgebra::*;
use crate::file_utils::read_vectors_from_file;

use hawklib::hawkkeygen::gen_f_g;
use hawklib::utils::rot_key;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use peak_alloc::PeakAlloc;

static PEAK_ALLOC: PeakAlloc = PeakAlloc;


pub fn run_hpp_attack(t: usize, n: usize) {
    // runs the HPP attack against Hawk
    //
    // Inputs: t digital signatures each signed with a secret key B
    // Outputs: Columns of B, or failure

    // In real life attack one would obtain the signatures as 
    // sig = enc(s1) where s1 = (h1 - w1) / 2
    // and one would need to recompute s0 on the signer/attacker side.
    // However for this implementation we simply generate entire signatures on signer end 
    // and directly return the entire signature w which we use to deduce information about B

    // get the correct key for later comparison
    let (b, binv) = get_secret_key(t, n);
    // correct G
    let cg = (&binv.transpose() * &binv).map(|x| x as f64);

    // STEP 0: Generate samples
    // We need to generate a lot of samples to run the attack on
    let samples = generate_samples(t, n);
    println!("Samples collected...");

    // STEP 1: estimate covariance matrix. This step requires a lot of samples,
    // so hopefully we can employ some sort of trick like the HPP against NTRU to reduce
    // number of signatures needed
    let g = estimate_covariance_matrix(&samples);
    println!("Covariance matrix estimated...");

    // STEP 2: conversion from hidden parallelepiped to hidden hypercube.  
    // in this step we need a covariance matrix estimation from step 1. The better
    // the estimation in step two, the better conversion estimation we can do here.
    // Given matrix G, we want to compute L s.t. L^t L = G^-1, and thereafter 
    // multiply our signature samples on the right with this L 
    // We use the Nalgebra crate for representations of matrices and for procedures such
    // as Cholesky decomposition
    let (u, linv) = hypercube_transformation(samples, cg);
    println!("Samples transformed...");
    
    // STEP 3: Gradient Descent:
    // The final step is to do gradient descent on our (converted) samples to minimize the 
    // fourth moment, and consequently reveal a row/column from +/- B
    if let Some(sol) = gradient_descent(u, 0.7){
        let res = (linv * sol).map(|x| x.round() as i16);
        eprintln!("{res}");
    }

}

fn get_secret_key(t: usize, degree: usize) -> (DMatrix<i32>, DMatrix<i32>) {
    // gets the secret key for t samples degree n
    // provided there is only one file

    // get the private key
    let (_, pkey) = read_vectors_from_file(&format!("{t}vectors_deg{degree}"))
        .expect(&format!("Could not find file with length {t} and degree {degree}"));

    // get the matrix form of b inverse
    let (b, binv) = to_mat(&pkey);

    (b.map(|x| x as i32), binv.map(|x| x as i32))
}

fn to_mat(privkey: &(Vec<u8>, Vec<i64>, Vec<i64>)) -> (DMatrix<i64>, DMatrix<i64>) {
    // given private key, reconstruct entire secret matrix B and B inverse
    let (fgseed, bigf, bigg) = privkey;
    let n = bigf.len();

    // reconstruct f and g
    let (f, g) = gen_f_g(fgseed, bigf.len());

    // create the matrix form of B and B inverse
    let b = rot_key(&f, &g, &bigf, &bigg);
    let binv = rot_key(
        &bigg,
        &g.clone().iter().map(|&x| -x).collect(),
        &bigf.clone().iter().map(|&x| -x).collect(),
        &f,
    );

    // convert them to Nalgebra matrices
    let flatb: Vec<i64> = b.into_iter().flatten().collect();
    let flatbinv: Vec<i64> = binv.into_iter().flatten().collect();

    let b = DMatrix::from_column_slice(2 * n, 2 * n, &flatb);
    let binv = DMatrix::from_column_slice(2 * n, 2 * n, &flatbinv);
    (b, binv)
}

fn generate_samples(t: usize, degree: usize) -> DMatrix<i32> {
    // returns samples of length t for Hawk degree n
    // will return data from file if file exists.
    // TODO: Otherwise, create the samples in place

    // read from precomputed file
    let (signatures, _) = read_vectors_from_file(&format!("{t}vectors_deg{degree}"))
        .expect(&format!("Could not find file with length {t} and degree {degree}"));
    // TODO create new samples if file does not exist

    // convert the vectors into a nalgebra matrix

    // get dimensions
    let rows = signatures.len();
    let cols = signatures[0].len();
    // flatten the signature samples
    let sig_flat: Vec<i16> = signatures.into_iter().flatten().collect();

    // construct nalgebra matrix from this
    let signature_matrix = DMatrix::from_row_slice(
        rows, cols, &sig_flat
        );

    // convert to i32
    let signature_matrix = signature_matrix.map(|x| x as i32);

    // return the nalgebra matrix
    signature_matrix
}

fn estimate_covariance_matrix(samples: &DMatrix<i32>) -> DMatrix<f64> {
    // estimate covariance matrix BtB given samples 

    // which scalar depends on Hawk degree and std.dev of distribution
    let sigma = match samples.ncols() / 2 {
        256 => 2.00205824, // experimentally, from samples vs 2 * 1.01 theorethically
        512 => 2.0 * 1.278,
        1024 => 2.0 * 1.299,
        _ => 0.0,
    };

    // store number of samples
    let nrows: f64 = samples.nrows() as f64;

    // convert the samples to f64 so we can do division on them later
    // let samples = samples.map(|x| x as i32);

    // compute yty and convert to f64 after computation is done
    let yty: DMatrix<f64> = (samples.transpose() * samples).map(|x| x as f64);
    // divide each entry by (sigma^2 * num_samples) to scale properly, and round the result
    let g: DMatrix<f64> = (yty/(sigma.powi(2)*nrows)).map(|x| x.round());

    // TODO check if this approach is faster
    // let g: DMatrix<f64> = (samples.transpose()/sigma.powi(2))*(&samples / nrows).map(|x| x.round());

    g
}

fn hypercube_transformation(samples: DMatrix<i32>, g: DMatrix<f64>) -> (DMatrix<f64>, DMatrix<f64>) {
    // given samples and estimate of covariance matrix, return transformed 
    // samples from hidden parallelepiped onto hidden hypercube for easier 
    // analysis later
    // Also returns the l inverse so we don't have to recompute it later

    // take inverse of G
    let ginv = g
        .try_inverse()
        .expect("Couldn't take inverse :(");

    // compute L = Cholesky decomposition of g inverse
    let l = Cholesky::new(ginv).expect("Couldn't do Cholesky decomposition of ginv :(");

    // compute inverse of L
    let linv = l
        .l()
        .clone()
        .try_inverse()
        .expect("Couldn't take inverse of l");

    // make a copy of samples converted to f64 to be able to multiply them with L
    let samples_f64: DMatrix<f64> = samples.map(|x| x as f64);
    let u = samples_f64.clone() * l.l();
    (u, linv)
}

fn gradient_descent(samples: DMatrix<f64>, delta: f64) -> Option<DVector<f64>> {
    // performs gradient descent on hypercube samples 
    
    let n = samples.ncols();
    let mut rng = StdRng::seed_from_u64(34632394);
    // 1: choose w uniformly from unit sphere of R^n
    let mut w = get_rand_w(n, &mut rng); 
    // 2: compute approx. gradient of nabla_mom_4
    loop {
        let g = grad_mom4(&w, &samples); 
        // 3: compute w_new = w-delta*g
        let mut w_new = &w - (delta*g)*&w; 
        // 4: normalize w_new
        w_new = &w_new/w_new.norm();
        // 5.1: if 4th moment of w_new is greater than 4th moment of w, we have "overshot" and return w
        if mom4(&w_new, &samples) >= mom4(&w, &samples) {
            return Some(w);
        }
        // 5.2: otherwise set w to be w_new and goto 2
        else {
            w = w_new;
        }
    }
    // if a solution cannot be found
    None
}

fn get_rand_w(n: usize, rng: &mut StdRng) -> DVector<f64> {
    // inputs:
    //  - n: number of samples to produce
    //  - dist_bound: the range to sample in
    //  - rng: a pre-seeded StdRng instance
    //
    //  outputs some randomly generated w on the unit circle

    use rand::distributions::{Distribution, Uniform};
    // define uniform distribution
    let dist = Uniform::from(-1.0..1.0);

    // initialize empty vector to store the samples
    let mut rnd_bytes: Vec<f64> = Vec::with_capacity(n);

    // sample n times
    for _ in 0..n {
        rnd_bytes.push(dist.sample(rng));
    }

    // load random number into DVector
    let mut w = DVector::from_vec(rnd_bytes);

    // normalize the w vector
    w /= w.norm();
    w
}

fn mom4(w: &DVector<f64>, samples: &DMatrix<f64>) -> f64 {
    // estimate 4th moment given samples and vector w
    // compute <u,w>^4
    let temp: DVector<f64> = (samples * w).map(|x| x.powi(4));
    // compute mean of above, and return result
    let res = temp.sum() / w.len() as f64;
    res
}

fn grad_mom4(w: &DVector<f64>, samples: &DMatrix<f64>) -> f64 {
    // estimate gradient of 4th moment given samples and vector w
    // compute 4(<u, w>^3 * u)
    let temp = (samples * w).map(|x| x.powi(3));
    let temp2 = &temp.transpose() * samples;
    let res = temp2.sum()*4.0/ w.len() as f64;
    res
} 
