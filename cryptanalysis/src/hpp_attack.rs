use crate::file_utils::read_vectors_from_file;
use crate::gradient_search::{gradient_descent, gradient_ascent};
use nalgebra::*;

use hawklib::hawkkeygen::gen_f_g;
use hawklib::ntru_solve::ntrusolve;
use hawklib::utils::rot_key;

use crate::test_candidate_vec::test_candidate_vec;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use std::io::{stdout, Write};

use peak_alloc::PeakAlloc;

static PEAK_ALLOC: PeakAlloc = PeakAlloc;
static TOLERANCE: f64 = 1e-10;

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
    // Samples is a tx2n matrix
    let samples = generate_samples(t, n);

    // STEP 1: estimate covariance matrix. This step is automatically given by public key Q

    // STEP 2: conversion from hidden parallelepiped to hidden hypercube.
    // Given Q, we do cholesky decomposition to get L s.t. Q = LLt
    // multiply our signature samples on the right with this L
    println!("U dim: {}, {}", samples.nrows(), samples.ncols());
    let (u, linv, c) = hypercube_transformation(samples, q, &binv);
    println!("U' dim: {}, {}", u.nrows(), u.ncols());
    // let total_num_elements = (u.nrows() * u.ncols()) as f64;
    // let mean = u.iter().sum::<f64>() / total_num_elements;
    //
    // let variance = u.iter()
    //     .map(|&x| {
    //         let diff = x-mean;
    //         diff*diff
    //     })
    // .sum::<f64>() / total_num_elements;
    //
    // let kurtosis = u.iter()
    //     .map(|&x| {
    //         let diff = (x-mean).powi(4);
    //         diff
    //     })
    // .sum::<f64>() / total_num_elements;
    //
    // println!("Mean: {}", mean);
    // println!("Var:  {}", variance);
    // println!("Kur:  {}", kurtosis);

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
        if let Some(sol) = gradient_descent(&u, &c) {
            res_gradient_descent = (&linv * &sol).map(|x| x.round() as i32);
                if vec_in_key(&res_gradient_descent, &binv) {
                    println!("FOUND! Result is in key based on direct checking");
                    break;
                } 
                measure_res(&res_gradient_descent, &binv);
                println!("Norm of res from descent: {}", res_gradient_descent.map(|x| x as f64).norm());
                println!("Norm of col0: {}", col0.map(|x| x as f64).norm());
                println!("Norm of coln: {}", coln.map(|x| x as f64).norm());
                println!("Result not in key... \n");
            }
    }

    println!("\nDoing gradient ascent...");
    retries = 0;
    while retries < 10{
        retries += 1;
        if let Some(sol) = gradient_ascent(&u, &c) {
            res_gradient_ascent = (&linv * &sol).map(|x| x.round() as i32);
                if vec_in_key(&res_gradient_ascent, &binv) {
                    println!("FOUND! Result is in key based on direct checking");
                    break;
                }
                measure_res(&res_gradient_ascent, &binv);
                println!("Norm of res from ascent: {}", res_gradient_ascent.map(|x| x as f64).norm());
                println!("Norm of col0: {}", col0.map(|x| x as f64).norm());
                println!("Norm of coln: {}", coln.map(|x| x as f64).norm());
                println!("Result not in key... \n");
            }
    }
}

fn hypercube_transformation(
    samples: DMatrix<i32>,
    q: DMatrix<f64>,
    skey: &DMatrix<i32>,
) -> (DMatrix<f64>, DMatrix<f64>, DMatrix<f64>) {
    // given samples and estimate of covariance matrix, return transformed
    // samples from hidden parallelepiped onto hidden hypercube for easier
    // analysis later
    // Also returns the l inverse so we don't have to recompute it later

    // get theoretical sigma here for scaling

    let sigma = match q.nrows() / 2{
        256 => 2.0 * 1.001, // 1.01 in theory but in practice it is measured to be 1.001
        512 => 2.0 * 1.278,
        _ =>   2.0 * 1.299,
    };

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
    let res = (&l.l().transpose() / sigma) * samples_f64;
    let c = &l.l().transpose() * skey.map(|x| x as f64);
    // eprintln!("{:.1}", &c * &c.transpose());
    // println!("Current mem usage: {} gb", PEAK_ALLOC.current_usage_as_gb());

    (res, linv, c)
}


fn vec_in_key(vec: &DVector<i32>, key: &DMatrix<i32>) -> bool {

    // Check if the vector exists as a column in the matrix
    let as_column = key.column_iter().any(|col| col == *vec);

    // Check if the negative vector exists as a column in the matrix
    let as_column_neg = key.column_iter().any(|col| -col == *vec);

    as_column || as_column_neg

}


fn measure_res(res: &DVector<i32>, binv: &DMatrix<i32>) {
    // given a solution, measure how far the solution is away from each column of the secret key
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    (0..res.len()).into_iter().for_each(|i| {
        let col: DVector<f64> = binv.column(i).into_owned().map(|x| x as f64);
        
        let diff = (col.abs() -res.map(|x| x as f64).abs()).norm();
        if diff < min { min = diff };
        if diff > max { max = diff };
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

    // println!("f: {:?}", f);
    // println!("g: {:?}", g);
    // println!("F: {:?}", bigf);
    // println!("G: {:?}", bigg);

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

