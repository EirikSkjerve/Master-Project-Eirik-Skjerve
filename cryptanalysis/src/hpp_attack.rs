use crate::file_utils::read_vectors_from_file;
use nalgebra::*;

use hawklib::hawkkeygen::gen_f_g;
use hawklib::ntru_solve::ntrusolve;
use hawklib::utils::rot_key;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use std::io::{stdout, Write};

use peak_alloc::PeakAlloc;

static PEAK_ALLOC: PeakAlloc = PeakAlloc;
static TOLERANCE: f64 = 1e-5;
static DELTA: f64 = 0.01;

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

    // get the public key Q
    let q = get_public_key(t, n).unwrap().map(|x| x as f64);

    println!("Running HPP attack with {t} samples against Hawk{n}");

    // STEP 0: Generate samples
    // We need to generate a lot of samples to run the attack on
    // samples is a tx2n matrix
    // TODO if samples does not exist they need to be generated
    let samples = generate_samples(t, n);
    // println!("Samples collected...");

    // println!("Current mem usage: {} gb", PEAK_ALLOC.current_usage_as_gb());

    // STEP 1: estimate covariance matrix. This step is automatically given by public key Q

    // STEP 2: conversion from hidden parallelepiped to hidden hypercube.
    // in this step we need a covariance matrix estimation from step 1. The better
    // the estimation in step two, the better conversion estimation we can do here.
    // Given matrix G, we want to compute L s.t. L^t L = G^-1, and thereafter
    // multiply our signature samples on the right with this L
    let (u, linv) = hypercube_transformation(samples, q, &binv);
    println!("Samples transformed...");

    // STEP 3: Gradient Descent:
    // The final step is to do gradient descent on our (converted) samples to minimize the
    // fourth moment, and consequently reveal a row/column from +/- B

    let col0 = binv.column(0);
    let coln = binv.column(n);

    let mut res_gradient_descent: DVector<i32> = DVector::zeros(2*n);
    let mut res_gradient_ascent: DVector<i32> = DVector::zeros(2*n);

    println!("\nDoing gradient descent...");
    let mut retries = 0;
    while retries < 10{
        retries += 1;
        if let Some(sol) = gradient_descent(&u, DELTA) {
            res_gradient_descent = (&linv * &sol).map(|x| x.round() as i32);
                if !vec_in_key(&res_gradient_descent, &binv) {
                    println!("Result not in key...");
                    continue;
                }
                println!("FOUND!!!");
                break;
            }
    }

    println!("\nDoing gradient ascent...");
    retries = 0;
    while retries < 10{
        retries += 1;
        if let Some(sol) = gradient_ascent(&u, DELTA) {
            res_gradient_ascent = (&linv * &sol).map(|x| x.round() as i32);
                if !vec_in_key(&res_gradient_ascent, &binv) {
                    println!("Result not in key...");
                    continue;
                }
                println!("FOUND!!!");
                break;
            }
    }
    println!("Norm of res from descent: {}", res_gradient_descent.map(|x| x as f64).norm());
    println!("Norm of res from ascent: {}", res_gradient_ascent.map(|x| x as f64).norm());
    println!("Norm of col0: {}", col0.map(|x| x as f64).norm());
    println!("Norm of coln: {}", coln.map(|x| x as f64).norm());
    measure_res(&res_gradient_descent, &binv);
    measure_res(&res_gradient_ascent, &binv);

    // TODO
    // for each res, we can try to construct matrix by using vec- and rot
    // functions, construct B' and see if B't B' is close to Q
    let (res1, res2) = (res_gradient_descent, res_gradient_ascent);
}

fn hypercube_transformation(
    samples: DMatrix<i32>,
    q: DMatrix<f64>,
    skey: &DMatrix<i32>,
) -> (DMatrix<f64>, DMatrix<f64>) {
    // given samples and estimate of covariance matrix, return transformed
    // samples from hidden parallelepiped onto hidden hypercube for easier
    // analysis later
    // Also returns the l inverse so we don't have to recompute it later

    // compute L = Cholesky decomposition of Q
    let l = Cholesky::new(q).expect("Couldn't do Cholesky decomposition of ginv");

    // compute inverse of Lt for later transformation back to parallelepiped
    let linv = l
        .l()
        .transpose()
        .clone()
        .try_inverse()
        .expect("Couldn't take inverse of l")
        .map(|x| x as f64);

    println!("Cholesky decomposition complete.");

    // println!("Current mem usage: {} gb", PEAK_ALLOC.current_usage_as_gb());
    let mut samples_f64: DMatrix<f64> = samples.map(|x| x as f64);
    // to guarantee that no extra memory is taken by this
    std::mem::drop(samples);
    // println!("Current mem usage: {} gb", PEAK_ALLOC.current_usage_as_gb());
    // multiply samples with L
    // samples = l.l().transpose() * &samples;
    let res = &l.l().transpose() * samples_f64;
    // println!("Current mem usage: {} gb", PEAK_ALLOC.current_usage_as_gb());

    (res, linv)
}

fn gradient_descent(samples: &DMatrix<f64>, delta: f64) -> Option<DVector<f64>> {
    // performs gradient descent on hypercube samples

    let n = samples.nrows();
    let mut rng = StdRng::seed_from_u64(rand::random::<u64>());
    // let mut rng = StdRng::seed_from_u64(1234);
    let mut num_iter = 0;
    // 1: choose w uniformly from unit sphere of R^n
    let mut w = get_rand_w(n, &mut rng);

    loop {
        num_iter += 1;
        println!("\nOn iteration {num_iter}...");
        // 2: compute approx. gradient of nabla_mom_4
        // println!("Computing gradient of 4th moment...");
        let g = grad_mom4(&w, &samples);

        // println!("Computing new w...");
        // 3: compute w_new = w-delta*g
        // eprintln!("g: {g}");
        let mut w_new = &w - (delta * g);

        // 4: normalize w_new
        w_new = &w_new / w_new.norm();

        // 5.1: if 4th moment of w_new is greater than 4th moment of w, we have "overshot" and return w
        // or if there is practically no change since last iteration
        let wnew_mom4 = mom4(&w_new, &samples);
        let w_mom4 = mom4(&w, &samples);
        println!("Mom4(w)    : {w_mom4}");
        println!("Mom4(w_new): {wnew_mom4}");
        if wnew_mom4 >= w_mom4 
            || (&w_new - &w).norm() <= 1e-3
        {
            println!("Returned in {num_iter} iterations!");
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

fn gradient_ascent(samples: &DMatrix<f64>, delta: f64) -> Option<DVector<f64>> {
    // performs gradient descent on hypercube samples

    let n = samples.nrows();
    let mut rng = StdRng::seed_from_u64(rand::random::<u64>());
    // let mut rng = StdRng::seed_from_u64(1234);

    let mut num_iter = 0;
    // 1: choose w uniformly from unit sphere of R^n
    let mut w = get_rand_w(n, &mut rng);

    loop {
        num_iter += 1;
        println!("\nOn iteration {num_iter}...");
        // 2: compute approx. gradient of nabla_mom_4
        // println!("Computing gradient of 4th moment...");
        let g = grad_mom4(&w, &samples);

        // println!("Computing new w...");
        // 3: compute w_new = w-delta*g
        let mut w_new = &w + (delta * g);

        // 4: normalize w_new
        w_new = &w_new / w_new.norm();
        // 5.1: if 4th moment of w_new is lower than 4th moment of w, we have "overshot" and return w
        // or if there is practically no change since last iteration
        // println!("{}", (&w_new - &w).abs().norm());
        let wnew_mom4 = mom4(&w_new, &samples);
        let w_mom4 = mom4(&w, &samples);
        println!("Mom4(w)    : {w_mom4}");
        println!("Mom4(w_new): {wnew_mom4}");
        if wnew_mom4 <=w_mom4 
            || (&w_new - &w).norm() <= 0.0001
        {
            println!("Returned in {num_iter} iterations!");
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

fn vec_in_key(vec: &DVector<i32>, key: &DMatrix<i32>) -> bool {
    // Check if the vector exists as a row in the matrix
    let as_row = key.row_iter().any(|row| row == vec.transpose());

    // Check if the vector exists as a column in the matrix
    let as_column = key.column_iter().any(|col| col == *vec);

    // Check if the negative vector exists as a row in the matrix
    let as_row_neg = key.row_iter().any(|row| -row == vec.transpose());

    // Check if the negative vector exists as a column in the matrix
    let as_column_neg = key.column_iter().any(|col| -col == *vec);

    as_row || as_column || as_row_neg || as_column_neg
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
    let dist = Uniform::from(-10.0..10.0);

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
    let temp: DVector<f64> = (samples.transpose() * w).map(|x| x.powi(4));
    // compute mean of above, and return result
    let res = temp.sum() / samples.ncols() as f64;
    res
}

fn grad_mom4(w: &DVector<f64>, samples: &DMatrix<f64>) -> DVector<f64> {
    let n = w.nrows();
    let t = samples.ncols();
    // estimate gradient of 4th moment given samples and vector w
    // compute 4(<u, w>^3 * u)
    // println!("Shape of U: ({}, {})", samples.nrows(), samples.ncols());
    // println!("Shape of w: ({}, {})", w.nrows(), w.ncols());
    // dot product
    // power of 3 to each entry
    // println!("Shape of uw: ({}, {})", uw3.nrows(), uw3.ncols());
    let uw3: DVector<f64> = (samples.transpose() * w).map(|x| x.powi(3));
    // let uw3u: DVector<f64> = (4.0 * (samples * uw3) / samples.nrows() as f64);
    // println!("Shape of uw3u: ({}, {})", uw3u.nrows(), uw3u.ncols());
    // let g: DVector<f64> = (4.0 * (samples * uw3) / samples.ncols() as f64);
    let uw3u = DMatrix::from_fn(n, t, |i, j| uw3[j] * samples[(i, j)]);
    let g = 4.0*uw3u.column_mean();
    // println!("Shape of g: ({}, {})", g.nrows(), g.ncols());
    g
}

fn measure_res(res: &DVector<i32>, binv: &DMatrix<i32>) {
    // given a solution, measure how far the solution is away from each column of the secret key
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    (0..res.len()).into_iter().for_each(|i| {
        let col: DVector<i32> = binv.column(i).into_owned();
        
        let dot = col.dot(res);
        let cos_theta = dot as f64 / (col.map(|x| x as f64).norm()*res.map(|x| x as f64).norm());
        let angle = cos_theta.acos().to_degrees();

        if angle > max { max = angle }
        if angle < min { min = angle }
    });

    println!("Min: {min} \nMax: {max}");
}

fn is_orthogonal(matrix: &DMatrix<f64>) -> bool {
    let identity = DMatrix::identity(matrix.ncols(), matrix.ncols());
    let qt_q = matrix.transpose() * matrix;
    let diff = (&qt_q - identity).norm();
    diff < TOLERANCE
}

fn is_orthonormal(matrix: &DMatrix<f64>) -> bool {
    is_orthogonal(matrix)
        && matrix
            .column_iter()
            .all(|col| (col.norm() - 1.0).abs() < TOLERANCE)
}

fn get_secret_key(t: usize, degree: usize) -> (DMatrix<i32>, DMatrix<i32>) {
    // gets the secret key for t samples degree n
    // provided there is only one file

    // get the private key
    let (_, pkey, _) = read_vectors_from_file(&format!("{t}vectors_deg{degree}")).expect(&format!(
        "Could not find file with length {t} and degree {degree}"
    ));
    // println!("F: {:?} \nG: {:?}", pkey.1, pkey.2);
    // get the matrix form of b inverse
    let (b, binv) = to_mat(&pkey);

    (b.map(|x| x as i32), binv.map(|x| x as i32))
}

fn get_public_key(t: usize, degree: usize) -> Option<DMatrix<i64>> {
    // get the public key
    let (_, _, (q00, q01)) = read_vectors_from_file(&format!("{t}vectors_deg{degree}")).expect(
        &format!("Could not find file with length {t} and degree {degree}"),
    );

    // use ntrusolve to get q01 and q11, and convert it into a DMatrix
    if let Some((q10, q11)) = ntrusolve(&q00, &q01) {
        let q = rot_key(&q00, &q10, &q01, &q11);
        let flatq: Vec<i64> = q.into_iter().flatten().collect();
        let q_mat = DMatrix::from_row_slice(2 * degree, 2 * degree, &flatq);
        return Some(q_mat);
    }
    None
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

    let b = DMatrix::from_row_slice(2 * n, 2 * n, &flatb);
    let binv = DMatrix::from_row_slice(2 * n, 2 * n, &flatbinv);
    (b, binv)
}

fn generate_samples(t: usize, degree: usize) -> DMatrix<i32> {
    // returns samples of length t for Hawk degree n
    // will return data from file if file exists.
    // TODO: Otherwise, create the samples in place

    // read from precomputed file
    let (signatures, _, _) = read_vectors_from_file(&format!("{t}vectors_deg{degree}")).expect(
        &format!("Could not find file with length {t} and degree {degree}"),
    );
    // TODO create new samples if file does not exist

    // convert the vectors into a nalgebra matrix

    // get dimensions
    let rows = signatures[0].len();
    let cols = signatures.len();

    // flatten the signature samples
    let sig_flat: Vec<i16> = signatures.into_iter().flatten().collect();

    // construct nalgebra matrix from this
    let signature_matrix = DMatrix::from_column_slice(rows, cols, &sig_flat);

    // convert to i32
    let signature_matrix = signature_matrix.map(|x| x as i32);

    // return the nalgebra matrix
    signature_matrix
}

