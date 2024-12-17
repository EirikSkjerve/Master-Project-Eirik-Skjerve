use hawklib::hawkkeygen::{gen_f_g, hawkkeygen};
use hawklib::hawksign::hawksign_total;
use hawklib::utils::rot_key;

use nalgebra::*;
use rand::Rng;

// type Matrix = DMatrix<i64>;

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

pub fn generate_t_signatures(t: usize, n: usize) {
    // create t signatures with hawk degree n

    // TODO check that rng seeds are big enough in hawksign w.r.t. number of samples

    assert!(n == 256 || n == 512 || n == 1024);

    // generate a keypair
    let (privkey, _) = hawkkeygen(n);

    // compute matrix version of secret key b and b inverse
    let (_, binv) = to_mat(&privkey);

    // compute g as b^-1t b^-1, covariance matrix of b inverse
    // let actual_g = binv.clone() * binv.transpose();
    let actual_g = binv.transpose() * binv;

    // create t messages
    let mut messages: Vec<Vec<u8>> = Vec::with_capacity(t);
    println!("Generating {t} messages...");
    for _ in 0..t {
        messages.push(get_random_bytes(100));
    }

    // create collection of t signatures corresponding to the above messages
    println!("Generating {t} signatures...");
    let mut signatures: Vec<Vec<i64>> = Vec::with_capacity(t);
    for i in 0..t {
        // sign each message
        // now each signature is on the form w = B^-1 * x
        signatures.push(hawksign_total(&privkey, &messages[i], n).0);
    }

    // convert the signatures into an t times 2n matrix
    // y contains our public samples distributed over P(B^-1)
    let sig_flat: Vec<i64> = signatures.into_iter().flatten().collect();
    let y = DMatrix::from_row_slice(t, 2 * n, &sig_flat);

    println!("Estimating covariance matrix...");
    let approx_g = estimate_cov_mat(&y);
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
    // let b = DMatrix::from_row_slice(2 * n, 2 * n, &flatb);
    // let binv = DMatrix::from_row_slice(2 * n, 2 * n, &flatbinv);

    let b = DMatrix::from_column_slice(2 * n, 2 * n, &flatb);
    let binv = DMatrix::from_column_slice(2 * n, 2 * n, &flatbinv);
    (b, binv)
}

fn estimate_cov_mat(y: &DMatrix<i64>) -> DMatrix<i64> {
    let sigma = match y.ncols() / 2 {
        256 => 2.0*1.01,
        512 => 2.0*1.278,
        1024 => 2.0*1.299,
        _ => 0.0,
    };

    let y_f = y.map(|x| x as f64);

    let g_approx = (1.0 / (sigma.powi(2))) * (y_f.transpose() * y_f) / y.nrows() as f64;

    let g_approx = g_approx.map(|x| x.round() as i64);
    g_approx
}

// gives a measure of the difference between two matrices
fn mat_dist(a_mat: &DMatrix<i64>, b_mat: &DMatrix<i64>) {
    let mut num_diff = 0;
    let mut sum_diff: i64 = 0;
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

    let avg_diff = sum_diff as f64 / num_diff as f64;
    println!("Matrices have different elements: {}", num_diff);
    println!("Average difference between elements: {}", avg_diff);
}
