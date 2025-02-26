use crate::file_utils::read_vectors_from_file;
use crate::gradient_search::{gradient_ascent, gradient_descent};
use crate::hpp_attack::measure_res;
use nalgebra::*;

use hawklib::hawkkeygen::gen_f_g;
use hawklib::utils::rot_key;

use hawklib::ntru_solve::ntrusolve;

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
static MAX_RETRIES: usize = 120;

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
    println!("Key generated");
    // eprintln!("B: {b}");
    // STEP 0: Generate samples
    // We need to generate a lot of samples to run the attack on
    // samples is a tx2n matrix

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
    (0..t).into_par_iter().for_each(|_|{
        let sig = hawk_sim_sign(n, &binv);
        signatures.lock().unwrap().push(hawk_sim_sign(n, &binv));
        pb.inc(1);
    });

    pb.finish_with_message("Signatures generated");

    let signatures = Arc::try_unwrap(signatures).expect("").into_inner().unwrap();
    let samples = DMatrix::from_columns(&signatures);
    let mut samples: DMatrix<f64> = samples.map(|x| x as f64);

    // STEP 1: estimate covariance matrix. We use Q=BtB for this

    let (linv, c) =
        hypercube_transformation(&mut samples, q.map(|x| x as f64), &binv.map(|x| x as i32));

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

    // run loop until number of max retries is set
    while retries < MAX_RETRIES {
        // increment counter
        retries += 1;

        // initialize empty result vector
        let mut res: Option<DVector<f64>> = None;

        // do both ascent and descent
        // Im desperado
        // if kurtosis is less than 3 we need to minimize
        if retries%2==0{
            println!("\nDoing gradient descent...");
            res = gradient_descent(&samples, correct_solution.as_ref());
        }

        // if kurtosis is greater than 3 we need to maximize
        if retries%2==1 {
            println!("\nDoing gradient ascent...");
            res = gradient_ascent(&samples, correct_solution.as_ref());
        }

        // multiply result vector with L inverse on the left to obtain solution as row in B
        // inverse
        let solution = (&linv * res.unwrap()).map(|x| x.round() as i32);

        // check directly if solution is in the actual secret key
        if vec_in_key(&solution, &binv.map(|x| x as i32)) {
            println!("FOUND! Result is in key based on direct checking");
            return;
        }


        // do a measurement of the result vector up against secret key if it was not the correct one
        measure_res(&solution, &binv.map(|x| x as i32));
        println!(
            "Norm of res from gradient search: {}",
            solution.map(|x| x as f64).norm()
        );
        println!("Norm of col0: {}", col0.map(|x| x as f64).norm());
        println!("Norm of coln: {}", coln.map(|x| x as f64).norm());
        println!("Result not in key... \n");
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

    // modify in place
    // this function is slow but avoids extra allocation

    // this method is fast but nalgebra matrix multiplication makes an extra allocation
    // for the matrices involved
    *samples = ((&l.l().transpose() / SIGMA) * &*samples);

    let c = &l.l().transpose() * skey.map(|x| x as f64);

    (linv, c)
}

fn vec_in_key(vec: &DVector<i32>, key: &DMatrix<i32>) -> bool {
    // Check if the vector exists as a column in the matrix
    let as_column = key.column_iter().any(|col| col == *vec);

    // Check if the negative vector exists as a column in the matrix
    let as_column_neg = key.column_iter().any(|col| -col == *vec);

    as_column || as_column_neg
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
