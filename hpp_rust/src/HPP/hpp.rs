use na::*;
use nalgebra as na;

use std::time::Instant;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use rand::distributions::{Distribution, Uniform};
// use rand_distr::{Normal, Distribution};

use crate::HPP::gradient_descent;

const NUM_SAMPLES: usize = 2000; 
const N: usize = 8;

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
    let dist = Uniform::from(-(dist_bound as f64)..(dist_bound as f64));
    // let normal = Normal::new(0.0, 1.0).unwrap();
    // initialize empty vector to store the samples
    let mut rnd_bytes: Vec<f64> = Vec::with_capacity(n);

    // sample n times and return the vector
    for _ in 0..n {
        rnd_bytes.push(dist.sample(rng));
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

fn check_v_approximation(v: &DMatrix<i32>, vapprox: &DMatrix<i32>) -> bool {
    // checks if one matrix is row permutation of rows of +/- rows of other matrix
    // return true/false based on this, and prints the permutation

    // displays the matrices (unsuitable for large parameters)
    // eprintln!("Approximation of rows of +- V: {vapprox}");
    // eprintln!("Actual V: {v}");

    // initialize an empty mapping vector
    let mut mapping: Vec<(usize, usize)> = Vec::with_capacity(N);

    // go through all rows in v
    for i in 0..N {
        let ra = v.row(i);

        // to through all rows in v~
        for j in 0..N {
            let rb = vapprox.row(j);

            // check for row equality a=b
            if ra == rb {
                // add position to map
                mapping.push((j, 0));
                break;
            }

            // check for row equality a=-b
            if ra == -rb {
                // add position to map
                mapping.push((j, 1));
                break;
            }
        }
    }

    println!("Matching rows: {}", mapping.len());

    // if mapping is "full", the matrices should be equal, because
    // rows in V are independent
    if mapping.len() == N {
        // print each mapping with correct sign
        // for (i, m) in mapping.iter().enumerate() {
        //     let (j, s) = m;
        //     if *s == 0 {
        //         println!("Row {} -> {}", i+1, j+1);
        //     }
        //     if *s == 1 {
        //         println!("Row {} -> -{}",i+1, j+1);
        //     }
        //
        // }
        println!("MATCH!");
        return true;
    }
    return false;
}

pub fn run_hpp_attack() {
    let start_pp = Instant::now();
    let entry_bound = 1;
    let dist_bound = 1;

    println!(
        "Running HPP attack on dimension {} with {} samples",
        N, NUM_SAMPLES
    );

    // generate some secret matrix V
    let sec_v_f = gen_sec_mat(N, entry_bound);

    let vtv = sec_v_f.transpose() * sec_v_f.clone();
    eprintln!("V^t V: {vtv}");

    // initialize
    // create random seed for rng
    let mut rng_seed = rand::thread_rng();
    let seed: u64 = rng_seed.gen();
    println!("Seed: {}", seed);

    let mut rng = StdRng::seed_from_u64(seed);

    // just hardcoding this for now
    let ex2 = 0.333;

    // empty vector storing samples
    let mut uni_samples = vec![];

    // generate a bunch of samples (that are uniformly distributed)
    // and multiply them with secret matrix v
    for i in 0..NUM_SAMPLES {
        let x = get_uni_slice_float(N, dist_bound, &mut rng);
        let x_vec = DVector::from_row_slice(&x);
        let y_vec = x_vec.transpose() * sec_v_f.clone();
        eprintln!("y_{i}: {y_vec}");
        uni_samples.push(y_vec);
    }

    // here the actual attack starts

    // now we have matrix Y
    let pub_y = DMatrix::from_rows(&uni_samples);

    // approximation of Gram Matrix
    let g_approx_f = (1.0 / ex2) * (pub_y.transpose() * pub_y.clone()) * (1.0 / NUM_SAMPLES as f64);

    // round the entries
    let g_approx = g_approx_f.map(|x| x.round());
    eprintln!("G: {g_approx}");

    // compute inverse of g
    let g_approx_inverse = g_approx
        .clone()
        .try_inverse()
        .expect("COULDN'T TAKE INVERSE OF G");

    // computing Cholesky decomposition of g⁻¹
    let l =
        Cholesky::new(g_approx_inverse).expect("COULDN'T COMPUTE CHOLESKY DECOMPOSITION OF Ginv");
    // extract lower triangular matrix l

    // get inverse of l
    // for some reason, taking l.inverse(), i.e. the nalgebra method for inverse directly on
    // the decomposition l would not work. Instead we retrieve the lowe triangular matrix L by .l()
    // and get the inverse from it like a normal square matrix
    let linv = l
        .l()
        .clone()
        .try_inverse()
        .expect("COULDN'T TAKE INVERSE OF L");

    let start = Instant::now();
    let u = pub_y * l.l(); // this should technically be divided by dist-bound
                           // pub_y * l / dist_bound as f64

    println!("Preprocessing used {:?}", start_pp.elapsed());

    let guess_sol = gradient_descent::gradient_descent(u, linv, 0.7);

    println!("gradient descent used: {:?}", start.elapsed());

    let sec_v = sec_v_f.map(|x| x.round() as i32);

    let rho_timer = Instant::now();
    check_v_approximation(&sec_v, &guess_sol);
    println!("Checking permutations took {:?}", rho_timer.elapsed());
}
