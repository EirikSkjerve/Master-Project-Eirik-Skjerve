use crate::file_utils::read_vectors_from_file;
use crate::gradient_search::{gradient_ascent, gradient_descent, gradient_descent_vanilla, gradient_ascent_vanilla};

use crate::collect_signatures::{collect_signatures_wh, get_random_bytes};

use crate::test_candidate_vec::test_candidate_vec;

use crate::compare_keys::compare_keys;

use nalgebra::*;

use hawklib::hawkkeygen::gen_f_g;
use hawklib::hawkkeygen::hawkkeygen;
use hawklib::hawksign::hawksign_total_h;
use hawklib::ntru_solve::ntrusolve;
use hawklib::utils::rot_key;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use rayon::prelude::*;
use indicatif::{ProgressBar, ProgressStyle};

use std::io::{stdout, Write};
use std::sync::{Arc, Mutex};
use std::time::{Instant, Duration};

use peak_alloc::PeakAlloc;

static PEAK_ALLOC: PeakAlloc = PeakAlloc;
static TOLERANCE: f64 = 1e-10;
static MAX_RETRIES: usize = 100;

static MU_Y: f64 = -0.25;
static VAR_Y: f64 = 1.1451;
static SIGMA_Y: f64 = 1.0701;
pub fn run_hpp_attack_yu(t: usize, n: usize) {

    println!("Running HPP attack with {t} samples against Hawk{n}");
    measure_d_distribution(t, n);
    let (mut samples, (b, binv), q) = generate_samples_and_keys(t, n).unwrap();

    // eprintln!("Before: {samples}");
    let (linv, c) = hypercube_transformation(&mut samples, q.map(|x| x as f64), &binv, &b);
    // eprintln!("After: {samples}");
    println!("Samples transformed");

    // do some measuring of moments
    let total_num_elements = (samples.nrows() * samples.ncols()) as f64;
    let mean = samples.iter().sum::<f64>() / total_num_elements;

    let variance = samples
        .iter()
        .map(|&x| {
            (x - mean).powi(2)
        })
        .sum::<f64>()
        / total_num_elements;

    let sigma = variance.sqrt();

    println!("Mean: {}", mean);
    println!("Var:  {}", variance);
    println!("Sigma: {}", sigma);

    // center the samples
    // samples = samples.map(|x| x - mean);

    // normalize the samples
    // samples /= sigma;


    let mean = samples.iter().sum::<f64>() / (total_num_elements * sigma);

    let variance = samples
        .iter()
        .map(|&x| {
            ((x/sigma) - mean).powi(2)
        })
        .sum::<f64>()
        / total_num_elements;

    let kurtosis = samples
        .iter()
        .map(|&x| {
            ((x / sigma) - mean).powi(4)
        })
        .sum::<f64>()
        / total_num_elements;

    println!("Mean: {}", mean);
    println!("Var:  {}", variance);
    // println!("Sigma: {}", sigma);
    println!("Kur:  {}", kurtosis);

    return;
    
    loop {
        let res = (&linv * gradient_ascent_vanilla(&samples).unwrap()).map(|x| x.round() as i32);

        // eprintln!("{res}");
        if vec_in_key(&res, &binv) {
            println!("FOUND! Result is in key based on direct checking");
            return;
        }
    }

}


fn hypercube_transformation(
    samples: &mut DMatrix<f64>,
    q: DMatrix<f64>,
    skey: &DMatrix<i32>,
    b: &DMatrix<i32>
) -> (DMatrix<f64>, DMatrix<f64>) {
    // given samples and and covariance matrix Q, return transformed
    // samples from hidden parallelepiped onto hidden hypercube
    // Also returns the l inverse so we don't have to recompute it later

    // get theoretical sigma here for scaling

    let start = Instant::now();

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
    *samples = (&l.l().transpose() * &*samples);
    // *samples = b.map(|x| x as f64) * &*samples;

    // for reference, compute the matrix C
    let c = &l.l().transpose() * skey.map(|x| x as f64);
    println!("Max usage so far: {} gb", PEAK_ALLOC.peak_usage_as_gb());

    println!("Times used for transforming samples: {:?}", start.elapsed());

    (linv, c)
}

pub fn measure_d_distribution(t: usize, n: usize) {

    let num_rounds = 10;
    for i in 0..num_rounds {
        let (privkey, _) = hawkkeygen(n, None);
        let pb = ProgressBar::new(3*t as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({per_sec})")
                .unwrap()
                .progress_chars("#>-"),
        );

        let mut mu = Arc::new(Mutex::new(0.0));

        (0..t).into_par_iter().for_each(|i| {
            let (_, _, d, _) = hawksign_total_h(&privkey, &get_random_bytes(100), n);
            let temp_mu: f64 = d.iter().map(|&x| (x as f64)).sum();
            *mu.lock().unwrap() += temp_mu;

            pb.inc(1);
        });

        let mu = *mu.lock().unwrap() / (2*t*n) as f64;

        let mut var = Arc::new(Mutex::new(0.0));
        (0..t).into_par_iter().for_each(|i| {
            let (_, _, d, _) = hawksign_total_h(&privkey, &get_random_bytes(100), n);
            let temp_var: f64 = d.iter().map(|&x| (x as f64 - mu).powi(2)).sum();
            *var.lock().unwrap() += temp_var;

            pb.inc(1);
        });

        let var = *var.lock().unwrap() / (2*t*n) as f64;
        let sigma = var.sqrt();

        let mut kur = Arc::new(Mutex::new(0.0));

        println!("Mu: {mu}");
        println!("Var: {var}");
        println!("Sigma: {sigma}");

        let mu = mu/sigma;

        (0..t).into_par_iter().for_each(|i| {
            let (_, _, d, _) = hawksign_total_h(&privkey, &get_random_bytes(100), n);
            let temp_kur: f64 = d.iter().map(|&x| ((x as f64 / sigma) - mu).powi(4)).sum();
            *kur.lock().unwrap() += temp_kur;

            pb.inc(1);

        });

        pb.finish_with_message("Estimation completed");

        let kur = *kur.lock().unwrap() / (2*t*n) as f64;


        println!("Norm. mu: {mu}");
        println!("Norm. kur: {kur}");
    }
}

fn vec_in_key(vec: &DVector<i32>, key: &DMatrix<i32>) -> bool {
    // Check if the vector exists as a column in the matrix
    let as_column = key.column_iter().any(|col| col == *vec);

    // Check if the negative vector exists as a column in the matrix
    let as_column_neg = key.column_iter().any(|col| -col == *vec);

    as_column || as_column_neg
}

pub fn measure_res(res: &DVector<i32>, binv: &DMatrix<i32>) -> (f64, f64){
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
    // println!("Min norm of diff: {min} \nMax norm of diff: {max}");
    (min, max)
}


pub fn to_mat(privkey: &(Vec<u8>, Vec<i64>, Vec<i64>)) -> (DMatrix<i64>, DMatrix<i64>) {
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

    let (signatures, privkey, pubkey) = collect_signatures_wh(t, degree);

    println!("{}, {}", signatures.len(), signatures[0].len());

    // flatten the ws samples
    let sig_flat: Vec<i16> = signatures.into_iter().flatten().collect();

    // construct nalgebra matrix from this
    let signature_matrix = DMatrix::from_column_slice(2*degree, t, &sig_flat);

    // convert to f64
    let signature_matrix = signature_matrix.map(|x| x as f64);

    println!("{}, {}", signature_matrix.nrows(), signature_matrix.ncols());

    // convert the private key (i.e. f and g) to the entire matrix B and B inverse
    let (b, binv) = to_mat(&privkey);
    // unpack 
    let (q00, q01) = pubkey;

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
