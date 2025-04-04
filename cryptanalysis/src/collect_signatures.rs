use hawklib::hawkkeygen::{gen_f_g, hawkkeygen};
use hawklib::hawksign::{hawksign_total, hawksign_total_h};
use hawklib::utils::rot_key;

use crate::file_utils::{read_vectors_from_file, write_vectors_to_file};

use nalgebra::*;
use rand::Rng;
use std::io::{stdout, Write};
use std::mem;
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;

use peak_alloc::PeakAlloc;
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

pub fn get_random_bytes(num_bytes: usize) -> Vec<u8> {
    // return num_bytes random bytes in a Vec<u8>
    //
    // example: get_random_bytes(4) -> vec![100,200,33,99]

    let mut res = Vec::with_capacity(num_bytes);
    let mut rng = rand::thread_rng();

    for _ in 0..num_bytes {
        res.push(rng.gen_range(0..255));
    }

    res
}

pub fn collect_signatures_wh(
    t: usize,
    n: usize,
) -> (
    Vec<Vec<i16>>,
    (Vec<u8>, Vec<i64>, Vec<i64>),
    (Vec<i64>, Vec<i64>),
) {
    // create t signatures with hawk degree n and store them in file

    println!("Collecting {t} signatures of Hawk degree {n}");

    let filename = format!("{t}vectors_deg{n}");

    // generate a keypair
    // let (privkey, pubkey) = hawkkeygen(n, Some(vec![1,3,3,7]));
    let (privkey, pubkey) = hawkkeygen(n, Some(vec![1, 2, 3]));

    let pb = ProgressBar::new(t as u64);

    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({per_sec})")
            .unwrap()
            .progress_chars("#>-"),
    );
    // create collection of t signatures corresponding to the above messages
    println!("Generating {t} signatures...");
    let mut ws: Arc<Mutex<Vec<Vec<i16>>>> = Arc::new(Mutex::new(Vec::with_capacity(t)));

    let mut zeros = Arc::new(Mutex::new(0));
    let mut ones = Arc::new(Mutex::new(0));
    (0..t).into_par_iter().for_each(|x| {
        let (_, h, _, w1) = hawksign_total_h(&privkey, &get_random_bytes(100), n);
        let w1: Vec<i16> = w1.iter().map(|&x| x as i16).collect();

        for i in h {
            if i == 0 {
                *zeros.lock().unwrap() += 1
            }
            if i == 1 {
                *ones.lock().unwrap() += 1
            }
        }

        ws.lock().unwrap().push(w1);

        pb.inc(1);
    });

    println!(
        "Zeros: {}\nOnes:  {}",
        *zeros.lock().unwrap(),
        *ones.lock().unwrap()
    );

    pb.finish_with_message("Completed");

    let ws = Arc::try_unwrap(ws)
        .expect("Could not unpack..")
        .into_inner()
        .unwrap();

    (ws, privkey, pubkey)
}

pub fn collect_signatures_par(
    t: usize,
    n: usize,
    write: bool,
) -> Option<(
    Vec<Vec<i16>>,
    (Vec<u8>, Vec<i64>, Vec<i64>),
    (Vec<i64>, Vec<i64>),
)> {
    let filename = format!("{t}vectors_deg{n}");

    // generate a keypair
    let (privkey, pubkey) = hawkkeygen(n, None);

    println!("Generating {t} signatures...");

    let mut signatures: Arc<Mutex<Vec<Vec<i16>>>> = Arc::new(Mutex::new(Vec::with_capacity(t)));
    // let mut stdout = Arc::new(Mutex::new(stdout()));
    let pb = ProgressBar::new(t as u64);

    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({per_sec})")
            .unwrap()
            .progress_chars("#>-"),
    );

    (0..t).into_par_iter().for_each(|i| {
        let sig: Vec<i16> = hawksign_total(&privkey, &get_random_bytes(100), n)
            .iter()
            .map(|&x| x as i16)
            .collect();

        signatures.lock().unwrap().push(sig);
        pb.inc(1);
    });

    pb.finish_with_message("Completed");

    let signatures_unpacked = Arc::try_unwrap(signatures)
        .expect("Could not unpack the signatures")
        .into_inner()
        .unwrap();

    if write {
        write_vectors_to_file(signatures_unpacked.to_vec(), privkey, pubkey, &filename);
        println!("\nWritten signatures to {}", filename);
        return None;
    } else {
        return Some((signatures_unpacked.to_vec(), privkey, pubkey));
    }
}

pub fn covariance_matrix_estimation(t: usize, n: usize) {
    // estimates covariance matrix of secret key
    // this function retrieves signatures from file

    assert!(n == 256 || n == 512 || n == 1024);

    println!("Estimating covariance matrix given {t} signature samples for Hawk degree {n}");

    // for nice printouts
    let mut stdout = stdout();

    // read from precomputed file
    let (signatures, privkey, pubkey) =
        read_vectors_from_file(&format!("{t}vectors_deg{n}")).expect("Could not read file");

    // compute matrix version of secret key b and b inverse
    let (_, binv) = to_mat(&privkey);

    // compute g as b^-1t b^-1, covariance matrix of b inverse
    let actual_g: DMatrix<i16> = (binv.transpose() * binv).map(|x| x as i16);

    println!("Current mem usage: {} MB", PEAK_ALLOC.current_usage_as_mb());

    // convert the signatures into an t times 2n matrix
    // y contains our public samples distributed over P(B^-1)
    // signatures variable is used up here and will not take any more memory
    let sig_flat: Vec<i16> = signatures.into_iter().flatten().collect();
    let y = DMatrix::from_row_slice(t, 2 * n, &sig_flat);
    // drop the flattened signatures for memory saving purposes
    std::mem::drop(sig_flat);

    println!(
        "Current mem usage after creating DMatrix: {} MB",
        PEAK_ALLOC.current_usage_as_mb()
    );

    println!("Estimating covariance matrix...");
    let approx_g = estimate_cov_mat(&y);
    println!("Estimated covariance matrix");
    mat_dist(&actual_g, &approx_g);
}

fn to_mat(privkey: &(Vec<u8>, Vec<i64>, Vec<i64>)) -> (DMatrix<i64>, DMatrix<i64>) {
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

fn estimate_cov_mat(y: &DMatrix<i16>) -> DMatrix<i16> {
    let sigma = match y.ncols() / 2 {
        256 => 2.00205824, // experimentally, from samples vs 2 * 1.01 theorethically
        512 => 2.0 * 1.278,
        1024 => 2.0 * 1.299,
        _ => 0.0,
    };

    let nrows: f32 = y.nrows() as f32;
    // first compute yt_y of i16 matrices to save as much memory as possible
    println!("max in y: {}", y.max());
    let y = y.map(|x| x as f32);
    println!(
        "Current mem usage after converting y to f32: {} MB",
        PEAK_ALLOC.current_usage_as_mb()
    );
    let g = (y.transpose() / sigma.powi(2)) * (&y / nrows);
    std::mem::drop(y);
    println!(
        "Current mem usage after computing yt_y: {} MB",
        PEAK_ALLOC.current_usage_as_mb()
    );
    println!("max in yt_y: {}", g.max());
    // now yt_y is nxn and will not take as much memory
    let g_r = g.map(|x| x.round() as i16);
    // return final result
    g_r
}

// gives a measure of the difference between two matrices
fn mat_dist(a_mat: &DMatrix<i16>, b_mat: &DMatrix<i16>) {
    let mut num_diff = 0;
    let mut sum_diff: i64 = 0;
    for i in 0..a_mat.nrows() {
        for j in 0..b_mat.nrows() {
            let a = a_mat[(i, j)];
            let b = b_mat[(i, j)];

            if a != b {
                // println!("{} != {}", a, b);
                num_diff += 1;
                sum_diff += (a as i64 - b as i64).abs();
            }
        }
    }

    let avg_diff = sum_diff as f64 / (a_mat.nrows() * a_mat.ncols()) as f64;
    println!("Matrices have different elements: {}", num_diff);
    println!("Average difference between elements: {}", avg_diff);
}
