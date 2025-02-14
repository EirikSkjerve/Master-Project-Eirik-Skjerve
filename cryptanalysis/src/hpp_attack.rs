use crate::file_utils::read_vectors_from_file;
use crate::gradient_search::{gradient_ascent, gradient_descent};

use crate::collect_signatures::collect_signatures_par;

use crate::test_candidate_vec::test_candidate_vec;

use nalgebra::*;

use hawklib::hawkkeygen::gen_f_g;
use hawklib::ntru_solve::ntrusolve;
use hawklib::utils::rot_key;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use rayon::prelude::*;

use std::io::{stdout, Write};
use std::sync::{Arc, Mutex};

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

    println!("Running HPP attack with {t} samples against Hawk{n}");
    // STEP 0: Generate samples
    // We need to generate a lot of samples to run the attack on
    // Samples is a tx2n matrix

    let (mut samples, (b, binv), q) = generate_samples_and_keys(t, n).unwrap();

    // STEP 1: estimate covariance matrix. This step is automatically given by public key Q

    // STEP 2: conversion from hidden parallelepiped to hidden hypercube.
    // Given Q, we do cholesky decomposition to get L s.t. Q = LLt
    // multiply our signature samples on the right with this L
    // to produce U which is distributed over L^t B^-1
    // now samples are modified in place so samples = U
    // the modified samples are modified as f32, but later cast into f64
    let (linv, c) = hypercube_transformation(&mut samples, q.map(|x| x as f64), &binv);
    println!("Samples transformed");

    let total_num_elements = (samples.nrows() * samples.ncols()) as f64;
    let mean = samples.iter().sum::<f64>() / total_num_elements;
    let variance = samples
        .iter()
        .map(|&x| {
            let diff = (x - mean).powi(2);
            diff
        })
        .sum::<f64>()
        / total_num_elements;

    let kurtosis = samples
        .iter()
        .map(|&x| {
            let diff = (x - mean).powi(4);
            diff
        })
        .sum::<f64>()
        / total_num_elements;

    println!("Mean: {}", mean);
    println!("Var:  {}", variance);
    println!("Kur:  {}", kurtosis);


    // STEP 3: Gradient Descent:
    // The final and main step is to do gradient descent on our (converted) samples to minimize the
    // fourth moment, and consequently reveal a row/column from +/- B

    let col0 = binv.column(0);
    let coln = binv.column(n);

    println!(
        "Current memory usage: {} mb",
        PEAK_ALLOC.current_usage_as_mb()
    );


    // random column for testing
    // let rand_index: usize = rand::thread_rng().gen_range(0..2*n);
    // let correct_solution: Option<DVector<f64>> = Some(c.column(5).into_owned());
    let correct_solution: Option<DVector<f64>> = None;

    // initialize retry counter
    let mut retries = 0;
    let max_retries = 10;
    // run loop until number of max retries is set
    while retries < max_retries {
        // increment counter
        retries += 1;

        // initialize empty result vector
        let mut res: Option<DVector<f64>> = None;

        // if kurtosis is less than 3 we need to minimize
        if kurtosis < 3.0 {
            println!("\nDoing gradient descent...");
            res = gradient_descent(&samples, correct_solution.as_ref());
        }

        // if kurtosis is greater than 3 we need to maximize
        if kurtosis > 3.0 {
            println!("\nDoing gradient ascent...");
            res = gradient_ascent(&samples, correct_solution.as_ref());
        }

        // check if result is a value
        if res.is_some() {

            // multiply result vector with L inverse on the left to obtain solution as row in B
            // inverse
            let solution = (&linv * res.clone().unwrap()).map(|x| x.round() as i32);

            // check directly if solution is in the actual secret key
            if vec_in_key(&solution, &binv) {
                println!("FOUND! Result is in key based on direct checking");
                return;
            }

            // do a measurement of the result vector up against secret key if it was not the correct one
            measure_res(&solution, &binv);
            println!(
                "Norm of res from gradient search: {}",
                solution.map(|x| x as f64).norm()
            );
            println!("Norm of col0: {}", col0.map(|x| x as f64).norm());
            println!("Norm of coln: {}", coln.map(|x| x as f64).norm());
            println!("Result not in key... \n");
        }
    }

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

    let sigma = match q.nrows() / 2 {
        256 => 2.0 * 1.001, // 1.01 in theory but in practice it is measured to be 1.001, unsure of
                           // which to choose 
        512 => 2.0 * 1.278,
        _ => 2.0 * 1.299,
    };

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

    // this method is fast but nalgebra matrix multiplication makes an extra allocation
    // for the matrices involved
    *samples = ((&l.l().transpose() / sigma) * &*samples);
    // *samples = &l.l().transpose() * &*samples;

    let c = &l.l().transpose() * skey.map(|x| x as f64);
    println!("Max usage so far: {} gb", PEAK_ALLOC.peak_usage_as_gb());

    (linv, c)
}


fn multiply_in_place_parallel(l: &DMatrix<f64>, m: &mut DMatrix<f64>) {
    let (rows_l, cols_l) = (l.nrows(), l.ncols());
    let (rows_m, cols_m) = (m.nrows(), m.ncols());

    assert_eq!(cols_l, rows_m, "Matrix dimensions do not match for multiplication!");

    // Wrap each row of `m` in a separate Mutex to allow parallel writes
    let row_mutexes: Vec<_> = m.row_iter_mut().map(|r| Arc::new(Mutex::new(r.into_owned()))).collect();

    (0..rows_l).into_par_iter().for_each(|i| {
        let mut row_result = vec![0.0; cols_m];  // Thread-local row buffer

        // Compute row `i` of the result
        for j in 0..cols_m {
            let mut sum = 0.0;
            for k in 0..cols_l {
                // Read safely from M (without locking the whole thing)
                let row_k = row_mutexes[k].lock().unwrap();
                sum += l[(i, k)] * row_k[j];
            }
            row_result[j] = sum;
        }

        // Write computed row back to `M`
        let mut row_lock = row_mutexes[i].lock().unwrap();
        row_lock.copy_from_slice(&row_result);
    });

    // Copy the modified row data back into `m`
    for (i, mut row) in m.row_iter_mut().enumerate() {
        let locked_row = row_mutexes[i].lock().unwrap();
        row.copy_from_slice(&locked_row.as_slice());
    }
}

//     let (rows_l, cols_l) = (lhs.nrows(), lhs.ncols());
//     let (rows_r, cols_r) = (rhs.nrows(), rhs.ncols());
//     // Temporary buffer to hold one row of the output at a time
//     let mut row_buffer = vec![0.0; rhs.ncols()];
//
//     // Iterate over rows of L
//     for i in 0..rows_l {
//         // Compute row `i` of the result and store in `row_buffer`
//         for j in 0..cols_r {
//             let mut sum = 0.0;
//             for k in 0..cols_l {
//                 sum += lhs[(i, k)] * rhs[(k, j)];
//             }
//             row_buffer[j] = sum;  // Store computed value
//         }
//
//         // Write computed row back to M (since we are modifying in place)
//         for j in 0..cols_r {
//             rhs[(i, j)] = row_buffer[j];
//         }
//     }
// }

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
    // eprintln!("{comb}");
    println!("Min norm of diff: {min} \nMax norm of diff: {max}");
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

fn generate_samples_and_keys(
    t: usize,
    degree: usize,
) -> Option<(DMatrix<f64>, (DMatrix<i32>, DMatrix<i32>), DMatrix<i64>)> {
    // check first if file exists
    let file_exists: bool = read_vectors_from_file(&format!("{t}vectors_deg{degree}")).is_ok();
    if file_exists {
        println!("{t}vectors_deg{degree}.bin exists");
    } else {
        println!("{t}vectors_deg{degree}.bin does not exist");
    }
    let (signatures, pkey, skey) = match file_exists {
        false => collect_signatures_par(t, degree, false).unwrap(),
        true => read_vectors_from_file(&format!("{t}vectors_deg{degree}")).unwrap(),
    };

    // get dimensions
    let rows = signatures[0].len();
    let cols = signatures.len();

    // flatten the signature samples
    let sig_flat: Vec<i16> = signatures.into_iter().flatten().collect();

    // construct nalgebra matrix from this
    let signature_matrix = DMatrix::from_column_slice(rows, cols, &sig_flat);

    // convert to f64
    let signature_matrix = signature_matrix.map(|x| x as f64);

    let (b, binv) = to_mat(&pkey);
    let (q00, q01) = skey;

    // use ntrusolve to get q01 and q11, and convert it into a DMatrix
    if let Some((q10, q11)) = ntrusolve(&q00, &q01) {
        let q = rot_key(&q00, &q10, &q01, &q11);
        let flatq: Vec<i64> = q.into_iter().flatten().collect();
        let q_mat = DMatrix::from_row_slice(2 * degree, 2 * degree, &flatq);

        return Some((
            signature_matrix,
            (b.map(|x| x as i32), binv.map(|x| x as i32)),
            q_mat,
        ));
    } else {
        return None;
    }
}

// fn gen_samples_and_keys_inner()

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
