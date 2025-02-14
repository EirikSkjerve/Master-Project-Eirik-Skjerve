use crate::file_utils::read_vectors_from_file;
use crate::gradient_search::{gradient_ascent, gradient_descent};
use nalgebra::*;

use hawklib::hawkkeygen::gen_f_g;
use hawklib::utils::rot_key;

use hawklib::ntru_solve::ntrusolve;

use rand::prelude::*;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use rand_distr::{Distribution, Normal, Uniform};

use std::io::{stdout, Write};

use peak_alloc::PeakAlloc;

static PEAK_ALLOC: PeakAlloc = PeakAlloc;
static TOLERANCE: f64 = 1e-12;
static SIGMA: f64 = 1.0;

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

    println!("Running HPP attack with {t} samples against Hawk{n} simulation");

    let (b, binv) = create_secret_matrix(n).unwrap();
    // eprintln!("B: {b} \nB^-1: {binv}");
    let b = b.map(|x| x as f64);
    let binv = binv.map(|x| x as f64);

    let q = &b.transpose() * &b;

    // STEP 0: Generate samples
    // We need to generate a lot of samples to run the attack on
    // samples is a tx2n matrix

    let mut samples_cols: Vec<DVector<f64>> = Vec::new();
    let mut rng = StdRng::seed_from_u64(rand::random::<u64>());
    for i in 0..t {
        let x = DVector::from_vec(get_norm_slice_float(2 * n, &mut rng)).map(|x| x.round());
        // let x = DVector::from_vec(get_uniform_slice_float(2*n, &mut rng));
        let y: DVector<f64> = &binv * x;
        samples_cols.push(y);
    }

    let mut samples = DMatrix::<f64>::from_columns(&samples_cols);

    // STEP 1: estimate covariance matrix. We use Q=BtB for this

    println!("U dim: {}, {}", samples.nrows(), samples.ncols());
    let (linv, c) =
        hypercube_transformation(&mut samples, q, &binv.map(|x| x as i32));
    println!("U' dim: {}, {}", samples.nrows(), samples.ncols());

    // // STEP 3: Gradient Descent:
    // // The final step is to do gradient descent on our (converted) samples to minimize the
    // // fourth moment, and consequently reveal a row/column from +/- B
    // let mut num_iter = 0;
    // loop {
    //     num_iter += 1;
    //     println!("Doing gradient descent...");
    //     if let Some(sol) = gradient_descent(&samples, &binv, None) {
    //         let res = (&linv * &sol).map(|x| x.round() as i32);
    //         if vec_in_key(&res, &binv.map(|x| x.round() as i32)) {
    //             println!("FOUND!");
    //             eprintln!("{res}");
    //             eprintln!("{binv}");
    //             println!("Total iterations: {num_iter}");
    //             break;
    //         }
    //         eprintln!("Res: {res}");
    //         println!("Norm of res: {}", res.map(|x| x as f64).norm());
    //         println!("Norm of col0: {}", binv.column(0).norm());
    //         println!("Norm of coln: {}", binv.column(n).norm());
    //         // eprintln!("{binv}");
    //     }
    //
    //     if num_iter == 10 {
    //         break;
    //     }
    // }
    //
    // num_iter = 0;
    // loop {
    //     println!("Doing gradient ascent...");
    //     if let Some(sol) = gradient_ascent(&samples, &binv, None) {
    //         let res = (&linv * &sol).map(|x| x.round() as i32);
    //         if vec_in_key(&res, &binv.map(|x| x.round() as i32)) {
    //             println!("FOUND!");
    //             eprintln!("{res}");
    //             eprintln!("{binv}");
    //             println!("Total iterations: {num_iter}");
    //             break;
    //         }
    //         eprintln!("Res: {res}");
    //         println!("Norm of res: {}", res.map(|x| x as f64).norm());
    //         println!("Norm of col0: {}", binv.column(0).norm());
    //         println!("Norm of coln: {}", binv.column(n).norm());
    //     }
    //     if num_iter == 10 {
    //         break;
    //     }
    // }
}

fn hypercube_transformation(
    samples: &mut DMatrix<f64>,
    q: DMatrix<f64>,
    skey: &DMatrix<i32>,
) -> (DMatrix<f64>, DMatrix<f64>) {
    // given samples and estimate of covariance matrix, return transformed
    // samples from hidden parallelepiped onto hidden hypercube for easier
    // analysis later
    // Also returns the l inverse so we don't have to recompute it later

    // get theoretical sigma here for scaling

    // compute L = Cholesky decomposition of Q
    let l = Cholesky::new(q).expect("Couldn't do Cholesky decomposition of ginv");

    // compute inverse of Lt for later transformation back to parallelepiped
    let linv = l
        .l()
        .transpose()
        .clone()
        .try_inverse()
        .expect("Couldn't take inverse of l");

    println!("Cholesky decomposition complete.");

    println!("Max usage so far: {} gb", PEAK_ALLOC.peak_usage_as_gb());
    // modify in place
    // this function is slow but avoids extra allocation

    // this method is fast but nalgebra matrix multiplication makes an extra allocation
    // for the matrices involved
    *samples = ((&l.l().transpose() / SIGMA) * &*samples);
    println!("Current mem usage: {} gb", PEAK_ALLOC.current_usage_as_gb());

    let c = &l.l().transpose() * skey.map(|x| x as f64);
    println!("Max usage so far: {} gb", PEAK_ALLOC.peak_usage_as_gb());
    println!("Current mem usage: {} gb", PEAK_ALLOC.current_usage_as_gb());

    (linv, c)
}
fn measure_res(res: &DVector<i32>, binv: &DMatrix<i32>) {
    // given a solution, measure how far the solution is away from each column of the secret key
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    let mut min_index = 0;
    (0..res.len()).into_iter().for_each(|i| {
        let col: DVector<f64> = binv.column(i).into_owned().map(|x| x as f64);

        let diff = (col.abs() - res.map(|x| x as f64).abs()).norm();
        if diff < min {
            min = diff;
            min_index = i
        };
        if diff > max {
            max = diff
        };
    });

    let comb = DMatrix::from_columns(&[res.column(0), binv.column(min_index)]);
    eprintln!("{comb}");
    println!("Min norm of diff: {min} \nMax norm of diff: {max}");
}
fn create_secret_matrix(n: usize) -> Option<(DMatrix<i64>, DMatrix<i64>)> {
    // here we try and replicate a HAWK matrix for smaller degree
    let mut i = 0;
    let mut rng = StdRng::seed_from_u64(rand::random());
    loop {
        i += 1;

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

fn create_secret_matrix_plain(n: usize) -> Option<(DMatrix<i64>, DMatrix<i64>)> {
    let mut j = 0;
    loop {
        let mut rng = StdRng::seed_from_u64(9876 + 13 * j);
        let mut columns: Vec<Vec<i64>> = Vec::new();
        for i in (0..2 * n) {
            columns.push(get_uniform_slice_fixed(2 * n, &mut rng));
        }
        let columnsflat: Vec<i64> = columns.into_iter().flatten().collect();
        let mat = DMatrix::from_row_slice(2 * n, 2 * n, &columnsflat);

        // TODO check for rank of mat and return if rank == 2n
        let svd = SVD::new(mat.map(|x| x as f64), true, true);
        let rank = svd
            .singular_values
            .iter()
            .filter(|&&sigma| sigma > TOLERANCE)
            .count();
        if rank == 2 * n {
            let mat_inverse = mat
                .map(|x| x as f64)
                .try_inverse()
                .unwrap()
                .map(|x| x.round() as i64);
            return Some((mat, mat_inverse));
        }
    }
}
// create plain uniform random secret matrix

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

pub fn get_uniform_slice_float(n: usize, rng: &mut StdRng) -> Vec<f64> {
    let dist = Uniform::new(-4.0 * SIGMA, 4.0 * SIGMA);

    let mut rnd_bytes: Vec<f64> = Vec::with_capacity(n);

    (0..n)
        .into_iter()
        .for_each(|_| rnd_bytes.push(dist.sample(rng).round()));
    rnd_bytes
}

pub fn get_uniform_slice_fixed(n: usize, rng: &mut StdRng) -> Vec<i64> {
    let dist = Uniform::new(-1, 1);

    let mut rnd_bytes: Vec<i64> = Vec::with_capacity(n);

    (0..n)
        .into_iter()
        .for_each(|_| rnd_bytes.push(dist.sample(rng)));
    rnd_bytes
}


fn vec_in_key(vec: &DVector<i32>, key: &DMatrix<i32>) -> bool {
    // Check if the vector exists as a column in the matrix
    let as_column = key.column_iter().any(|col| col == *vec);

    // Check if the negative vector exists as a column in the matrix
    let as_column_neg = key.column_iter().any(|col| -col == *vec);

    as_column || as_column_neg
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
