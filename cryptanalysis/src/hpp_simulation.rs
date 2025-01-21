use crate::file_utils::read_vectors_from_file;
use nalgebra::*;

use hawklib::hawkkeygen::gen_f_g;
use hawklib::utils::rot_key;

use hawklib::ntru_solve::ntrusolve;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use rand_distr::{Distribution, Normal, Uniform};

use std::io::{stdout, Write};

use peak_alloc::PeakAlloc;

static PEAK_ALLOC: PeakAlloc = PeakAlloc;
static TOLERANCE: f64 = 1e-10;
static DELTA: f64 = 0.1;
static SIGMA: f64 = 2.0;

pub fn run_hpp_sim(t: usize, n: usize) {
    // runs the HPP attack against Hawk
    //
    // Inputs: t digital signatures each signed with a secret key B
    // Outputs: Columns of B, or failure

    // In real life attack one would obtain the signatures as
    // sig = enc(s1) where s1 = (h1 - w1) / 2
    // and one would need to recompute s0 on the signer/attacker side.
    // However for this implementation we simply generate entire signatures on signer end
    // and directly return the entire signature w which we use to deduce information about B

    let (b, binv) = create_secret_matrix(n).unwrap();
    let b = b.map(|x| x as f64);
    let binv = binv.map(|x| x as f64);
    // eprintln!("binv: {b}");

    // Perform Singular Value Decomposition
    let svd = SVD::new(binv.clone(), true, true);

    // Compute the rank by counting singular values above the tolerance
    let rank = svd
        .singular_values
        .iter()
        .filter(|&&x| x > TOLERANCE)
        .count();

    println!("Rank of B: {}", rank);

    println!("Running HPP attack with {t} samples against Hawk{n}");

    let cg = &binv * &binv.transpose();
    // eprintln!("{cg}");

    // STEP 0: Generate samples
    // We need to generate a lot of samples to run the attack on
    // samples is a tx2n matrix

    let mut samples_rows: Vec<RowDVector<f64>> = Vec::new();
    for i in 0..t {
        let mut rng = StdRng::seed_from_u64(1003 + i as u64);
        let x = DVector::from_vec(get_norm_slice_float(2 * n, &mut rng));
        let y: DVector<f64> = &binv * x;
        samples_rows.push(y.transpose());
    }

    let samples = DMatrix::<f64>::from_rows(&samples_rows);
    println!(
        "Samples is a {} by {} matrix.",
        samples.nrows(),
        samples.ncols()
    );

    // STEP 1: estimate covariance matrix. This step requires a lot of samples,
    // so hopefully we can employ some sort of trick like the HPP against NTRU to reduce
    // number of signatures needed

    let g = estimate_covariance_matrix(&samples.map(|x| x as i32));
    mat_dist(&cg, &g);
    // eprintln!("{g}");
    println!("Covariance matrix estimated...");

    // STEP 2: conversion from hidden parallelepiped to hidden hypercube.
    // in this step we need a covariance matrix estimation from step 1. The better
    // the estimation in step two, the better conversion estimation we can do here.
    // Given matrix G, we want to compute L s.t. L^t L = G^-1, and thereafter
    // multiply our signature samples on the right with this L
    // We use the Nalgebra crate for representations of matrices and for procedures such
    // as Cholesky decomposition
    // CURRENTLY ONLY USING THE CORRECT G INSTEAD OF ESTIMATE

    let (u, linv) = hypercube_transformation(samples, cg, &binv);
    println!("Samples transformed...");

    // // STEP 3: Gradient Descent:
    // // The final step is to do gradient descent on our (converted) samples to minimize the
    // // fourth moment, and consequently reveal a row/column from +/- B
    // println!("Doing gradient descent...");
    // if let Some(sol) = gradient_descent(&u, DELTA) {
    //     let res = (&linv * &sol).map(|x| x.round() as i32);
    //     eprintln!("Result: {res}");
    //     println!("Is res in key? \n{} \n", vec_in_key(&res, &binv));
    // }
}

fn create_secret_matrix(n: usize) -> Option<(DMatrix<i64>, DMatrix<i64>)> {
    // here we try and replicate a HAWK matrix for smaller degree
    let mut i = 0;
    loop {
        let mut rng = StdRng::seed_from_u64(16789 + i);

        // generate uniformly distributed vectors f and g
        let (f, g) = (
            get_norm_slice_float(n, &mut rng),
            get_norm_slice_float(n, &mut rng),
        );

        let f_i64: Vec<i64> = f.iter().map(|x| x.round() as i64).collect();
        let g_i64: Vec<i64> = g.iter().map(|x| x.round() as i64).collect();

        if let Some((bigf, bigg)) = ntrusolve(&f_i64, &g_i64) {
            return Some(to_mat((f_i64, g_i64, bigf, bigg)));
        }
    }
}

fn to_mat(privkey: (Vec<i64>, Vec<i64>, Vec<i64>, Vec<i64>)) -> (DMatrix<i64>, DMatrix<i64>) {
    // given private key, construct entire secret matrix B and B inverse

    let (f, g, bigf, bigg) = privkey;
    let n = f.len();

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

/// returns n f64 floats uniformly distributed on -dist_bound..dist_bound
/// requires a pre-seeded StdRng instance
pub fn get_norm_slice_float(n: usize, rng: &mut StdRng) -> Vec<f64> {
    // inputs:
    //  - n: number of samples to produce
    //  - dist_bound: the range to sample in
    //  - rng: a pre-seeded StdRng instance

    // define upper and lower bound for the sampling
    let dist = Normal::new(0.0, SIGMA).unwrap();
    // initialize empty vector to store the samples
    let mut rnd_bytes: Vec<f64> = Vec::with_capacity(n);

    // sample n times and return the vector
    for _ in 0..n {
        rnd_bytes.push(dist.sample(rng));
    }
    rnd_bytes
}

fn estimate_covariance_matrix(samples: &DMatrix<i32>) -> DMatrix<f64> {
    // estimate covariance matrix BtB given samples

    // which scalar depends on Hawk degree and std.dev of distribution
    let sigma = match samples.ncols() / 2 {
        256 => 2.00205824, // experimentally, from samples vs 2 * 1.01 theorethically
        512 => 2.0 * 1.278,
        1024 => 2.0 * 1.299,
        _ => SIGMA,
    };

    // store number of samples
    let nrows: f64 = samples.nrows() as f64;

    // compute yty and convert to f64 after computation is done
    let yty: DMatrix<f64> = (samples.transpose() * samples).map(|x| x as f64);
    // divide each entry by (sigma^2 * num_samples) to scale properly, and round the result
    println!("{sigma}");
    let g: DMatrix<f64> = (yty).map(|x| (x / (sigma.powi(2) * nrows)).round());

    let is_positive_definite = g.clone().cholesky().is_some();
    println!("Is G positive definite? {is_positive_definite}");

    // TODO check if this approach is faster
    // let g: DMatrix<f64> = (samples.transpose()/sigma.powi(2))*(&samples / nrows).map(|x| x.round());

    // Perform Singular Value Decomposition
    let svd = SVD::new(g.clone(), true, true);

    // Compute the rank by counting singular values above the tolerance
    let rank = svd
        .singular_values
        .iter()
        .filter(|&&x| x > TOLERANCE)
        .count();

    eprintln!("G: {g}");
    println!("Rank of G: {}", rank);
    assert_eq!(g, g.transpose());

    g
}

fn hypercube_transformation(
    samples: DMatrix<f64>,
    g: DMatrix<f64>,
    skey: &DMatrix<f64>,
) -> (DMatrix<f64>, DMatrix<f64>) {
    // given samples and estimate of covariance matrix, return transformed
    // samples from hidden parallelepiped onto hidden hypercube for easier
    // analysis later
    // Also returns the l inverse so we don't have to recompute it later

    // take inverse of G
    let ginv = g.try_inverse().expect("Couldn't take inverse :(");
    // compute L = Cholesky decomposition of g inverse
    let l = Cholesky::new(ginv).expect("Couldn't do Cholesky decomposition of ginv :(");

    // compute inverse of L
    let linv = l
        .l()
        .clone()
        .try_inverse()
        .expect("Couldn't take inverse of l");

    // make a copy of samples converted to f64 to be able to multiply them with L
    let samples_transformed: DMatrix<f64> = samples * l.l();

    let c = skey.transpose() * l.l();
    // let cct = &c * &c.transpose();
    // eprintln!("{cct:.2}");
    if is_orthogonal(&c) {
        println!("C is orthogonal!");
    } else {
        println!("C is not orthogonal...");
    }

    (samples_transformed, linv)
}

fn gradient_descent(samples: &DMatrix<f64>, delta: f64) -> Option<DVector<f64>> {
    // performs gradient descent on hypercube samples

    let n = samples.ncols();
    let mut rng = StdRng::seed_from_u64(332394);
    let mut num_iter = 0;
    // 1: choose w uniformly from unit sphere of R^n
    let mut w = get_rand_w(n, &mut rng);

    // 2: compute approx. gradient of nabla_mom_4
    loop {
        num_iter += 1;
        let g = grad_mom4(&w, &samples);
        // 3: compute w_new = w-delta*g
        let mut w_new = &w - (delta * g);
        // 4: normalize w_new
        w_new = &w_new / w_new.norm();
        // 5.1: if 4th moment of w_new is greater than 4th moment of w, we have "overshot" and return w
        if mom4(&w_new, &samples) >= mom4(&w, &samples) {
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
    // performs gradient descent on hypercube by maximizing 4th moment

    let n = samples.ncols();
    let mut rng = StdRng::seed_from_u64(34872114);
    let mut num_iter = 0;
    let mut stdout = stdout();
    // 1: choose w uniformly from unit sphere of R^n
    let mut w = get_rand_w(n, &mut rng);
    // 2: compute approx. gradient of nabla_mom_4
    let mut prevnorm = 0.0;
    loop {
        // print!("\rIterations: {num_iter}");
        // std::io::Write::flush(&mut std::io::stdout()).unwrap();
        num_iter += 1;
        let g = grad_mom4(&w, &samples);
        let curnorm = g.norm();
        // println!("\n|g|={}", curnorm);
        // eprintln!("g: {g:.2}");
        // 3: compute w_new = w-delta*g
        let mut w_new = &w + (delta * g);
        // 4: normalize w_new
        w_new = &w_new / w_new.norm();

        // println!("mom4(w_new)={}", mom4(&w_new, &samples));
        // println!("mom4(w)={}", mom4(&w, &samples));

        // if (curnorm-prevnorm).abs() < 1.0 {
        //     return Some(w)
        // }

        // 5.1: if 4th moment of w_new is smaller than 4th moment of w, we have "overshot" and return w
        if mom4(&w_new, &samples) <= mom4(&w, &samples) {
            std::io::Write::flush(&mut std::io::stdout()).unwrap();
            println!("\nReturned in {num_iter} iterations!");
            return Some(w);
        }
        // 5.2: otherwise set w to be w_new and goto 2
        else {
            w = w_new;
            prevnorm = curnorm;
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

    as_row || as_column
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

fn grad_mom4(w: &DVector<f64>, samples: &DMatrix<f64>) -> DVector<f64> {
    // estimate gradient of 4th moment given samples and vector w
    // compute 4(<u, w>^3 * u)

    // dot product
    let uw3: DVector<f64> = (samples * w).map(|x| x.powi(3));
    // power of 3 to each entry
    let uw3u: DVector<f64> =
        (4.0 * (uw3.transpose() * samples) / samples.nrows() as f64).transpose();
    uw3u
}

fn is_orthogonal(matrix: &DMatrix<f64>) -> bool {
    let identity = DMatrix::identity(matrix.ncols(), matrix.ncols());
    let qt_q = matrix.transpose() * matrix;
    let diff = (&qt_q - identity).norm();
    // eprintln!("{qt_q:.1}");
    println!("Diff: {diff}");
    diff < TOLERANCE
}

fn is_orthonormal(matrix: &DMatrix<f64>) -> bool {
    is_orthogonal(matrix)
        && matrix
            .column_iter()
            .all(|col| (col.norm() - 1.0).abs() < TOLERANCE)
}

// gives a measure of the difference between two matrices
fn mat_dist(a_mat: &DMatrix<f64>, b_mat: &DMatrix<f64>) {
    let mut num_diff = 0;
    let mut sum_diff: f64 = 0.0;
    for i in 0..a_mat.nrows() {
        for j in 0..b_mat.nrows() {
            let a = a_mat[(i, j)];
            let b = b_mat[(i, j)];

            if a != b {
                // println!("{} != {}", a, b);
                num_diff += 1;
                sum_diff += (a - b).abs();
            }
        }
    }

    let avg_diff = sum_diff as f64 / (a_mat.nrows() * a_mat.ncols()) as f64;
    println!("Matrices have different elements: {}", num_diff);
    println!("Average difference between elements: {}", avg_diff);
}
