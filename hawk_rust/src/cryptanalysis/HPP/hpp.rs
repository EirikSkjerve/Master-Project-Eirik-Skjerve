use na::*;
use nalgebra as na;
use std::collections::HashSet;
use std::time::Instant;

use rand::distributions::{Distribution, Uniform};
use rand::rngs::StdRng;
use rand::SeedableRng;

use crate::cryptanalysis::HPP::gradient_descent;

const NUM_SAMPLES: usize = 1000;
const N: usize = 5;

/// returns n i8 integers uniformly distributed on -entry_bound..entry_bound
pub fn get_uni_slice_int(n: usize, entry_bound: usize, seed: usize) -> Vec<i8> {
    // inputs:
    //  - n: number of samples to produce
    //  - entry_bound: the range to sample in
    //  - seed: a seed for the StdRng instance

    // define upper and lower bound for the sampling
    let bound = Uniform::from(-(entry_bound as i8)..(entry_bound as i8) + 1);
    // seed an rng
    let mut rng = StdRng::seed_from_u64(seed as u64);
    // initialize empty vector to store the samples
    let mut rnd_bytes: Vec<i8> = Vec::with_capacity(n);

    // sample n times and return the vector
    for _ in 0..n {
        rnd_bytes.push(bound.sample(&mut rng));
    }
    rnd_bytes
}

/// returns n f64 floats uniformly distributed on -dist_bound..dist_bound
/// requires a pre-seeded StdRng instance
pub fn get_uni_slice_float(n: usize, dist_bound: usize, rng: &mut StdRng) -> Vec<f64> {
    // inputs:
    //  - n: number of samples to produce
    //  - dist_bound: the range to sample in
    //  - rng: a pre-seeded StdRng instance

    // define upper and lower bound for the sampling
    let bound = Uniform::from(-(dist_bound as f64)..(dist_bound as f64));
    // initialize empty vector to store the samples
    let mut rnd_bytes: Vec<f64> = Vec::with_capacity(n);

    // sample n times and return the vector
    for _ in 0..n {
        rnd_bytes.push(bound.sample(rng));
    }
    rnd_bytes
}

/// generate a secret full rank, square [degree x degree] matrix V
/// with entries uniformly distributed on -entry_bound..entry_bound
fn gen_sec_mat(
    degree: usize,
    entry_bound: usize,
) -> Matrix<f64, Dyn, Dyn, VecStorage<f64, Dyn, Dyn>> {
    loop {
        // generate slice of uniformly random integers between -entry_bound..entry_bound
        let uni_slice = get_uni_slice_int(degree * degree, entry_bound, 42);
        // convert the integers to floats for later calculation
        let uni_slice_f: Vec<f64> = uni_slice.iter().map(|&x| x as f64).collect();

        // create a matrix from this slice
        let sec_v = DMatrix::from_column_slice(degree, degree, &uni_slice_f);

        // perform singular value decomposition to calculate rank of matrix
        let svd_v = sec_v.clone().svd(true, true);
        let rank = svd_v.rank(1e-10);

        // rerun the loop if the rank is not max
        if rank < degree {
            continue;
        }

        // return secret matrix v
        return sec_v;
    }
}

fn map_rows(
    vapprox: Vec<Matrix<i32, Const<1>, Dyn, VecStorage<i32, Const<1>, Dyn>>>,
    v: Matrix<i32, Dyn, Dyn, VecStorage<i32, Dyn, Dyn>>,
) {
    let mut v_copy = v.clone();
    let mut ctr = 0;
    for vp in vapprox.iter() {
        for i in 0..v_copy.nrows() {
            if vp == &v_copy.row(i) {
                ctr += 1;
                v_copy = v_copy.clone().remove_rows_at(&[i]);
                break;
            }
        }
    }
}

pub fn run_hpp_attack() {
    // let start = Instant::now();
    let entry_bound = 1;
    let dist_bound = 1;

    // using a fixed matrix V
    // generate some secret matrix V
    let sec_v_f = gen_sec_mat(N, entry_bound);
    let v_data: [i32; 5 * 5] = [
        1, 1, 0, 0, -1, -1, 1, 0, 1, -1, 0, 1, 0, -1, -1, 1, -1, 1, 1, 1, -1, -1, 0, 0, 0,
    ];

    let v_data_f: Vec<f64> = v_data.iter().map(|&x| x as f64).collect();
    let sec_v_f = DMatrix::from_row_slice(5, 5, &v_data_f);

    // initialize
    let mut rng = StdRng::seed_from_u64(99999);
    let mut ctr = 0;

    ctr += 1;

    // just hardcoding this for now
    let ex2 = 0.333;

    // empty vector storing samples
    let mut uni_samples = vec![];

    // generate a bunch of samples (that are uniformly distributed)
    // and multiply them with secret matrix v
    for _ in 0..NUM_SAMPLES {
        let x = get_uni_slice_float(N, dist_bound, &mut rng);
        let x_vec = DVector::from_row_slice(&x);
        let y_vec = x_vec.transpose() * sec_v_f.clone();
        uni_samples.push(y_vec);
    }

    // now we have matrix Y
    let pub_y = DMatrix::from_rows(&uni_samples);

    // approximation of Gram Matrix
    let g_approx_f = (1.0 / ex2) * (pub_y.transpose() * pub_y.clone()) * (1.0 / NUM_SAMPLES as f64);

    let g_approx = g_approx_f.map(|x| x.round());

    eprintln!("g: {g_approx}");
    let mut g_approx_inverse = g_approx.clone();
    match g_approx_inverse.clone().try_inverse() {
        Some(g_inv) => {
            g_approx_inverse = g_inv;
        }
        None => {
            println!("Not invertible");
        }
    }

    // computing Cholesky decomposition of g⁻¹
    let l = Cholesky::new(g_approx_inverse).expect("TEST");
    let l = l.l();

    // get inverse of l
    let linv = l.clone().try_inverse().expect("COULDN'T TAKE INVERSE");

    let start = Instant::now();
    let u = pub_y * l; // this should technically be divided by dist-bound
                       // pub_y * l / dist_bound as f64
    let guess_sol = gradient_descent::gradient_descent(u, linv, 0.75);

    println!("gradient descent used: {:?}", start.elapsed());

    eprintln!("{}", sec_v_f);
}
