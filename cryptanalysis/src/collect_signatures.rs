use hawklib::hawkkeygen::{gen_f_g, hawkkeygen};
use hawklib::hawksign::hawksign_total;
use hawklib::utils::rot_key;

use crate::file_utils::{write_vectors_to_file, read_vectors_from_file};

use nalgebra::*;
use rand::Rng;
use std::io::{stdout, Write};
use std::time::{Duration, Instant};
use std::mem;

use peak_alloc::PeakAlloc;
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

fn get_random_bytes(num_bytes: usize) -> Vec<u8> {
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

pub fn collect_signatures(t: usize, n: usize) {
    // create t signatures with hawk degree n and store them in file
    
    let filename = format!("{t}vectors_deg{n}");

    let mut stdout = stdout();

    // generate a keypair
    let (privkey, _) = hawkkeygen(n);

    // create t messages
    let mut messages: Vec<Vec<u8>> = Vec::with_capacity(t);
    println!("Generating {t} messages...");
    for _ in 0..t {
        messages.push(get_random_bytes(100));
    }

    // create collection of t signatures corresponding to the above messages
    println!("Generating {t} signatures...");
    let mut signatures: Vec<Vec<i16>> = Vec::with_capacity(t);
    for i in 0..t {
        // sign each message
        // now each signature is on the form w = B^-1 * x
        // convert to Vec<i16> to save a lot of memory
        let sig: Vec<i16> = hawksign_total(&privkey, &messages[i], n)
            .0
            .iter()
            .map(|&x| x as i16)
            .collect();
        signatures.push(sig);

        // Calculate and display progress
        if i % (t / 100) == 0 || i == t - 1 {
            let progress = (i as f64 / t as f64) * 100.0;
            print!("\rProgress: {:.0}%", progress);
            std::io::Write::flush(&mut std::io::stdout()).unwrap();
        }
    }

    write_vectors_to_file(signatures, privkey, &filename);
    println!("\nWritten signatures to {}", filename);

}

pub fn covariance_matrix_estimation(t: usize, n: usize) {
    // estimates covariance matrix of secret key
    // this function retrieves signatures from file

    assert!(n == 256 || n == 512 || n == 1024);

    println!("Estimating covariance matrix given {t} signature samples for Hawk degree {n}");

    // for nice printouts
    let mut stdout = stdout();

    // read from precomputed file
    let (signatures, privkey) = read_vectors_from_file(&format!("{t}vectors_deg{n}"))
        .expect("Could not read file");

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

    println!("Current mem usage after creating DMatrix: {} MB", PEAK_ALLOC.current_usage_as_mb());

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
    println!("Current mem usage after converting y to f32: {} MB", PEAK_ALLOC.current_usage_as_mb());
    let g = (1.0 / sigma.powi(2))*y.transpose()*&y / nrows;
    std::mem::drop(y);
    println!("Current mem usage after computing yt_y: {} MB", PEAK_ALLOC.current_usage_as_mb());
    println!("max in yt_y: {}", g.max());
    // now yt_y is nxn and will not take as much memory
    // let g_f = (1.0 / sigma.powi(2)) * yt_y.map(|x| x as f64 ) / nrows;
    // println!("max in g approx: {}", g_f.max());
    // round the result and convert back to i16
    // let g_r = g_f.map(|x| x.round() as i16);
    // return final result
    let g_r = g.map(|x| x.round() as i16);
    g_r

    // let y_f = y.map(|x| x as f32);
    //
    // // we don't need y in memory any more
    // std::mem::drop(y);
    // println!("Current mem usage after converting y to floats: {} MB", PEAK_ALLOC.current_usage_as_mb());
    //
    // let g_approx = (1.0 / sigma.powi(2)) * y_f.transpose() * y_f / nrows;
    // println!("max in g approx: {}", g_approx.max());
    // let g_approx = g_approx.map(|x| x.round() as i16);
    // g_approx
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
