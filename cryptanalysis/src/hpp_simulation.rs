use crate::file_utils::read_vectors_from_file;
use crate::gradient_search::{
    gradient_ascent, gradient_ascent_vanilla, gradient_descent, gradient_descent_vanilla,
};
use crate::hpp_attack::measure_res;
use nalgebra::*;

use hawklib::hawkkeygen::gen_f_g;
use hawklib::utils::rot_key;

use hawklib::ntru_solve::ntrusolve;

use crate::compare_keys::compare_keys_simulated;
use crate::hawk_sim::*;

use rand::prelude::*;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use rand_distr::{Distribution, Normal, Uniform};

use std::io::{stdout, Write};
use std::sync::{Arc, Mutex};

use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;

use peak_alloc::PeakAlloc;

static PEAK_ALLOC: PeakAlloc = PeakAlloc;
static TOLERANCE: f64 = 1e-12;
static SIGMA: f64 = 2.02;
static MAX_RETRIES: usize = 1000;

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
    let ((b, binv), q) = hawk_sim_keygen(n);
    let bt = b.transpose();
    println!("Key generated");
    // eprintln!("B: {b}");
    // STEP 0: Generate samples
    // We need to generate a lot of samples to run the attack on
    // samples is a tx2n matrix

    // compare_keys_simulated(n);
    // return;

    // let mut stdout = Arc::new(Mutex::new(stdout()));
    let pb = ProgressBar::new(t as u64);

    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({per_sec})")
            .unwrap()
            .progress_chars("#>-"),
    );
    let mut samples_cols: Vec<DVector<i64>> = Vec::new();
    let mut signatures: Arc<Mutex<Vec<DVector<i64>>>> = Arc::new(Mutex::new(Vec::with_capacity(t)));
    // let mut rng = StdRng::seed_from_u64(rand::random::<u64>());
    (0..t).into_par_iter().for_each(|_| {
        let sig = hawk_sim_sign(n, &binv);
        signatures.lock().unwrap().push(sig);
        pb.inc(1);
    });

    pb.finish_with_message("Signatures generated");

    let signatures = Arc::try_unwrap(signatures).expect("").into_inner().unwrap();
    let samples = DMatrix::from_columns(&signatures);
    // let samp_min = samples.iter().min().unwrap();
    // let samp_max = samples.iter().max().unwrap();
    // println!("Min: {samp_min}");
    // println!("Max: {samp_max}");
    // println!("Clipping values");
    // let mut samples: DMatrix<f64> = rectangular_clipping(&samples, SIGMA, n as f64).map(|x| x as f64);
    let mut samples: DMatrix<f64> = samples.map(|x| x as f64);

    // STEP 1: estimate covariance matrix. We use Q=BtB for this

    let (linv, c) =
        hypercube_transformation(&mut samples, q.map(|x| x as f64), &bt.map(|x| x as i32));

    // do some measuring of moments
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

    let col0 = binv.column(0);
    let coln = binv.column(n);

    println!(
        "Current memory usage: {} mb",
        PEAK_ALLOC.current_usage_as_mb()
    );

    // random column for testing
    // can be randomized which column index is chosen
    // can either be None or some column
    // let rand_index: usize = rand::thread_rng().gen_range(0..2*n);
    // let correct_solution: Option<DVector<f64>> = Some(c.column(5).into_owned());
    let correct_solution: Option<DVector<f64>> = None;

    // initialize retry counter
    let mut retries = 0;

    let mut avg_min = 0.0;
    let mut avg_max = 0.0;
    let mut tot_min = f64::INFINITY;
    let mut tot_max = f64::NEG_INFINITY;
    // run loop until number of max retries is set
    while retries < MAX_RETRIES {
        // increment counter
        retries += 1;
        println!("Iteration {retries}");

        // initialize empty result vector
        let mut res: Option<DVector<f64>> = None;

        let res = gradient_descent_vanilla(&samples);

        // if (kurtosis - 3.0) >= 0.0 {
        //     println!("Doing gradient ascent");
        //     res = gradient_ascent_vanilla(&samples, (kurtosis-3.0));
        //     // res = gradient_ascent(&samples, None);
        // }
        // if (kurtosis - 3.0) < 0.0 {
        //     println!("Doing gradient descent");
        //     res = gradient_descent_vanilla(&samples, (kurtosis-3.0));
        //     // res = gradient_descent(&samples, None);
        // }

        // multiply result vector with L inverse on the left to obtain solution as row in B
        // inverse
        let solution = (&linv * res.unwrap()).map(|x| x.round() as i32);

        // check directly if solution is in the actual secret key
        if vec_in_key(&solution, &bt.map(|x| x as i32)) {
            println!("FOUND! Result is in key based on direct checking");
            return;
        }

        // do a measurement of the result vector up against secret key if it was not the correct one
        let (min, max) = measure_res(&solution, &bt.map(|x| x as i32));
        avg_min += min / MAX_RETRIES as f64;
        avg_max += max / MAX_RETRIES as f64;
        if min < tot_min {
            tot_min = min
        }
        if max > tot_max {
            tot_max = max
        }
        // println!(
        //     "Norm of res from gradient search: {}",
        //     solution.map(|x| x as f64).norm()
        // );
        // println!("Norm of col0: {}", col0.map(|x| x as f64).norm());
        // println!("Norm of coln: {}", coln.map(|x| x as f64).norm());
        println!("Result not in key... \n");
    }
    // println!("Avg min: {avg_min} \n Avg max: {avg_max}");
    // println!("Total min: {tot_min} \nTotal max: {tot_max}");
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

    // compute L = Cholesky decomposition of Q^-1
    let l = Cholesky::new(q.clone().try_inverse().expect("Could not take inverse"))
        .expect("Couldn't do Cholesky decomposition of qinv");

    // compute inverse of Lt for later transformation back to parallelepiped
    let linv = l
        .l()
        .transpose()
        .clone()
        .try_inverse()
        .expect("Couldn't take inverse of l");

    println!("Cholesky decomposition complete.");

    // modify in place
    // this function is slow but avoids extra allocation

    // this method is fast but nalgebra matrix multiplication makes an extra allocation
    // for the matrices involved

    *samples = q * &*samples;
    // now samples are on the form w = B^T x
    // *samples = ((&l.l().transpose()) * &*samples) / SIGMA;
    *samples = ((&l.l().transpose()) * &*samples) / 10.0;
    // let min_value = samples.iter().cloned().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
    // let max_value = samples.iter().cloned().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
    // println!("Min: {min_value}");
    // println!("Max: {max_value}");

    let threshold = 6.0;
    *samples = rectangular_clipping(&samples, threshold);

    // println!("Before: {}", samples.ncols());
    //
    // // *samples = &*samples / threshold;
    // //
    // // // remove outlier samples to create a more cubic shape of the samples
    // *samples = norm_shaving(samples, threshold) / max_value;
    // println!("After: {}", samples.ncols());

    // *samples = ((&l.l().transpose()) * &*samples) / 10.0;

    let c = &l.l().transpose() * skey.map(|x| x as f64);
    // let ctc = &c.transpose()*&c;
    // let cct = &c*&c.transpose();
    // eprintln!("{ctc:.2}");
    // eprintln!("{cct:.2}");

    (linv, c)
}

fn norm_shaving(samples: &DMatrix<f64>, threshold: f64) -> DMatrix<f64> {
    let mut result = Vec::new();

    for c in 0..samples.ncols() {
        let col = samples.column(c); // Get the c-th column
        if col.norm() <= threshold {
            result.push(col.clone()); // If the norm is less than or equal to threshold, keep the column
        }
    }

    // Convert the result into a new matrix
    DMatrix::from_columns(&result)
}
/// Clips points to enforce a hypercubic shape in n-dimensional space.
/// Removes any column vector where any coordinate exceeds k * std.
///
/// # Arguments
/// * `w` - A DMatrix<i64> where each column is a point in n-dimensional space.
/// * 'threshold' - A float determining the upper and lower limit for which a vector should be
/// retained
///
/// # Returns
/// A DMatrix<i64> containing only the retained column vectors.
fn rectangular_clipping(w: &DMatrix<f64>, threshold: f64) -> DMatrix<f64> {
    // Collect columns that satisfy the threshold condition
    let retained_columns: Vec<_> = (0..w.ncols())
        .filter(|&j| {
            w.column(j)
                .iter()
                .all(|&val| (val.abs() as f64) < threshold)
        })
        .flat_map(|j| w.column(j).iter().copied().collect::<Vec<f64>>())
        .collect();

    // Convert retained columns back to a DMatrix<i64>
    if retained_columns.is_empty() {
        return DMatrix::<f64>::zeros(w.nrows(), 0); // Return empty matrix if nothing remains
    }

    let num_rows = w.nrows();
    let num_cols = retained_columns.len() / num_rows;
    DMatrix::from_vec(num_rows, num_cols, retained_columns)
}

/// Clips points based on Z-scores.
fn z_score_clipping(w: &DMatrix<f64>, k: f64) -> DMatrix<f64> {
    let mut retained_columns = Vec::new();

    for j in 0..w.ncols() {
        let mut valid = true;

        // Check each dimension of the column (each point)
        for i in 0..w.nrows() {
            let z_score = (w[(i, j)] - 0.0) / 2.0 * SIGMA;

            if z_score.abs() > k {
                valid = false;
                break;
            }
        }

        // If the point is valid (within the Z-score threshold), keep it
        if valid {
            let column = w.column(j).clone_owned();
            retained_columns.push(column);
        }
    }

    // Construct a new matrix with the valid columns
    let num_rows = w.nrows();
    let num_cols = retained_columns.len();
    DMatrix::from_columns(&retained_columns)
}

fn vec_in_key(vec: &DVector<i32>, key: &DMatrix<i32>) -> bool {
    // Check if the vector exists as a column in the matrix
    let as_column = key.column_iter().any(|col| col == *vec);

    // Check if the negative vector exists as a column in the matrix
    let as_column_neg = key.column_iter().any(|col| -col == *vec);

    as_column || as_column_neg
}

fn mat_dist(a_mat: &DMatrix<f64>, b_mat: &DMatrix<f64>) {
    // gives a measure of the difference between two matrices
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
