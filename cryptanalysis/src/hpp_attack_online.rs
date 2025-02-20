use hawklib::hawkkeygen::{gen_f_g, hawkkeygen};
use hawklib::hawksign::hawksign_total;
use hawklib::utils::rot_key;
use hawklib::ntru_solve::ntrusolve;

use nalgebra::*;

use rand_distr::{Distribution, Normal, Uniform};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use std::io::{stdout, Write};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use peak_alloc::PeakAlloc;

use rayon::prelude::*;

static PEAK_ALLOC: PeakAlloc = PeakAlloc;
static DELTA: f64 = 0.75;


pub fn hpp_attack_online(t: usize, n: usize) {

    // STEP 1: generate hawk keypair
    // compute L from public Q via cholesky decomposition

    let sigma = match n{
        256 => 2.0 * 1.001,
        512 => 2.0 * 1.278,
        _ => 2.0 * 1.299,
    };

    let keypair = hawkkeygen(n);
    let (privkey, pubkey) = keypair;
    let (b, binv) = to_mat_priv(&privkey);
    let q = to_mat_pub(&pubkey);
    // compute L = Cholesky decomposition of Q
    let l = Cholesky::new(q).expect("Couldn't do Cholesky decomposition of ginv");

    // compute inverse of Lt for later transformation back to parallelepiped
    let linv = l
        .l()
        .transpose()
        .clone()
        .try_inverse()
        .expect("Couldn't take inverse of l");

    let lt = &l.l().transpose() / sigma;

    let c = &l.l().transpose() * binv.map(|x| x as f64);

    // STEP 2: gradient search
    let candidate = (&linv * gradient_search("descent", &lt, &privkey, t)).map(|x| x.round() as i32);
    if vec_in_key(&candidate, &binv) {
        println!("FOUND!");
    }
}

fn gradient_search(
    dir: &str, 
    lt: &DMatrix<f64>, 
    privkey: &(Vec<u8>, Vec<i64>, Vec<i64>), 
    t: usize
    ) -> DVector<f64> {

    println!("Doing gradient {dir}");
    let n = privkey.1.len();

    let seed = rand::random::<u64>();

    let mut rng = StdRng::seed_from_u64(seed); 
    // random starting point
    let mut w = get_rand_w(2*n, &mut rng);    

    let mut counter = 0;
    // let mut mom4 = estimate_mom4_par(lt, t, privkey, &w);
    loop {
        println!("At iteration {}...", counter + 1);
        counter += 1;
        let (mut mom4, gradmom4) = estimate_mom4_gradmom4_par(lt, t, privkey, &w);
        let mut wnew = match dir {
            "ascent" => &w + DELTA*&gradmom4,
            "descent" => &w - DELTA*&gradmom4,
            _ => DVector::zeros(n),
        };

        let gradmom4_norm = &gradmom4.norm();
        println!("Norm of gradient g: {}", gradmom4_norm);
        // normalize wnew to put it back on the unit sphere
        wnew = wnew.normalize();

        // compute 4th moment of new w
        let mom4_new = estimate_mom4_par(&lt, t, privkey, &wnew);

        println!("Mom4(w)    : {mom4}");
        println!("Mom4(w_new): {mom4_new}");
        println!("Diff       : {}", mom4 - mom4_new);

        if (dir=="ascent" && mom4_new < mom4) || (dir=="descent" && mom4_new > mom4) {
            println!("Direction has changed. Possibly at extremum");
            return w;
        }

        w = wnew;
        mom4 = mom4_new;
    }

    DVector::<f64>::zeros(n)
    
}

fn estimate_mom4_gradmom4_par(
    lt: &DMatrix<f64>, 
    t: usize,
    privkey: &(Vec<u8>, Vec<i64>, Vec<i64>),
    w: &DVector<f64>
    ) -> (f64, DVector<f64>){

    let n = lt.ncols();
    let mut gradmom4 = Arc::new(Mutex::new(DVector::<f64>::zeros(n)));
    let mut mom4 = Arc::new(Mutex::new(0.0));

    // for each iteration, generate a signature and update mom4 and gradmom4 iteratively
    for i in 0..t {
        // generate random message
        let msg = get_random_bytes(30);
        // compute a Hawk signature given private key
        let signature: DVector<f64> = DVector::from_vec(hawksign_total(&privkey, &msg, n/2)).map(|x| x as f64);
        // transform signature by multiplying on the left with L^t
        let transig = lt * signature;
        
        // update mom4 and gradmom4
        *gradmom4.lock().unwrap() += transig.dot(w).powi(3) * &transig;
        *mom4.lock().unwrap() += transig.dot(w).powi(4);
    }


    let mut gradmom4 = gradmom4.lock().unwrap().clone();

    let mom4 = *mom4.lock().unwrap() / t as f64;

    // taking mean
    gradmom4 *= 4.0;
    gradmom4 /= t as f64;

    // return values
    (mom4, gradmom4)

}
fn estimate_gradmom4_par(
    lt: &DMatrix<f64>, 
    t: usize,
    privkey: &(Vec<u8>, Vec<i64>, Vec<i64>),
    w: &DVector<f64>
    ) -> DVector<f64>{

    let n = lt.ncols();
    let mut gradmom4 = Arc::new(Mutex::new(DVector::<f64>::zeros(n)));

    // for each iteration, generate a signature and update mom4 and gradmom4 iteratively
    for i in 0..t {
        // generate random message
        let msg = get_random_bytes(30);
        // compute a Hawk signature given private key
        let signature: DVector<f64> = DVector::from_vec(hawksign_total(&privkey, &msg, n/2)).map(|x| x as f64);
        // transform signature by multiplying on the left with L^t
        let transig = lt * signature;
        
        // update mom4 and gradmom4
        *gradmom4.lock().unwrap() += transig.dot(w).powi(3) * transig;
    }

    let mut gradmom4 = gradmom4.lock().unwrap().clone();

    // taking mean
    gradmom4 *= 4.0;
    gradmom4 /= t as f64;

    // return value
    gradmom4

}

fn estimate_gradmom4(
    lt: &DMatrix<f64>, 
    t: usize,
    privkey: &(Vec<u8>, Vec<i64>, Vec<i64>),
    w: &DVector<f64>
    ) -> DVector<f64>{

    let n = lt.ncols();
    let mut mom4: f64 = 0.0;
    let mut gradmom4: DVector<f64> = DVector::<f64>::zeros(n);

    // for each iteration, generate a signature and update mom4 and gradmom4 iteratively
    for i in 0..t {
        // generate random message
        let msg = get_random_bytes(30);
        // compute a Hawk signature given private key
        let signature: DVector<f64> = DVector::from_vec(hawksign_total(&privkey, &msg, n/2)).map(|x| x as f64);
        // transform signature by multiplying on the left with L^t
        let transig = lt * signature;
        
        // update mom4 and gradmom4
        gradmom4 += transig.dot(w).powi(3) * transig;
    }

    // taking mean
    gradmom4 *= 4.0;
    gradmom4 /= t as f64;

    // return value
    gradmom4

}

fn estimate_mom4_par(
    lt: &DMatrix<f64>, 
    t: usize,
    privkey: &(Vec<u8>, Vec<i64>, Vec<i64>),
    w: &DVector<f64>
    ) -> f64 {

    let n = lt.ncols();
    let mut mom4 = Arc::new(Mutex::new(0.0));

    (0..t).into_par_iter().for_each(|_| {
        // generate random message
        let msg = get_random_bytes(30);
        // compute a Hawk signature given private key
        let signature: DVector<f64> = DVector::from_vec(hawksign_total(&privkey, &msg, n/2)).map(|x| x as f64);
        // transform signature by multiplying on the left with L^t
        let transig = lt * signature;

        // update mom4 and gradmom4
        *mom4.lock().unwrap() += transig.dot(w).powi(4);
    });

    let mom4 = *mom4.lock().unwrap() / t as f64;
    mom4
}

fn estimate_mom4(
    lt: &DMatrix<f64>, 
    t: usize,
    privkey: &(Vec<u8>, Vec<i64>, Vec<i64>),
    w: &DVector<f64>
    ) -> f64 {

    let n = lt.ncols();
    let mut mom4: f64 = 0.0;

    // for each iteration, generate a signature and update mom4 and gradmom4 iteratively
    for i in 0..t {
        // generate random message
        let msg = get_random_bytes(30);
        // compute a Hawk signature given private key
        let signature: DVector<f64> = DVector::from_vec(hawksign_total(&privkey, &msg, n/2)).map(|x| x as f64);
        // transform signature by multiplying on the left with L^t
        let transig = lt * signature;

        // update mom4 and gradmom4
        mom4 += transig.dot(w).powi(4);
    }

    // taking mean of number
    mom4 /= t as f64;
    // return value
    mom4

}

fn get_rand_w(n: usize, rng: &mut StdRng) -> DVector<f64> {
    //  outputs some randomly generated w on the unit circle

    let sigma = match n / 2 {
        256 => 2.0 * 1.001,
        512 => 2.0 * 1.278,
        _ => 2.0 * 1.299,
    };
    // define uniform distribution
    let dist1 = Normal::new(0.0, 5.0).unwrap();
    let dist2 = Normal::new(0.0, sigma).unwrap();
    let dist3 = Uniform::from(-10.0..10.0);

    // initialize empty vector to store the samples
    let mut rnd_bytes: Vec<f64> = Vec::with_capacity(n);

    // sample n times
    // n/2 for each distribution
    for _ in 0..n/2 {
        rnd_bytes.push(dist3.sample(rng));
    }
    for _ in 0..n/2 {
        rnd_bytes.push(dist3.sample(rng));
    }

    // load random number into DVector
    let mut w = DVector::from_vec(rnd_bytes);
    // normalize the w vector
    w = w.normalize();
    w
}
fn to_mat_pub(pubkey: &(Vec<i64>, Vec<i64>)) -> DMatrix<f64> {

    // use ntrusolve to get q01 and q11, and convert it into a DMatrix
    let (q00, q01) = pubkey;
    let n = q00.len();
    let (q10, q11) = ntrusolve(&q00, &q01).unwrap();
    let q_nest = rot_key(&q00, &q10, &q01, &q11);
    let flatq: Vec<i64> = q_nest.into_iter().flatten().collect();
    let q = DMatrix::from_row_slice(2 * n, 2 * n, &flatq).map(|x| x as f64);
    q
}

fn to_mat_priv(privkey: &(Vec<u8>, Vec<i64>, Vec<i64>)) -> (DMatrix<i32>, DMatrix<i32>) {
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

    let b = DMatrix::from_row_slice(2 * n, 2 * n, &flatb).map(|x| x as i32);
    let binv = DMatrix::from_row_slice(2 * n, 2 * n, &flatbinv).map(|x| x as i32);
    (b, binv)
}

fn vec_in_key(vec: &DVector<i32>, key: &DMatrix<i32>) -> bool {
    // Check if the vector exists as a column in the matrix
    let as_column = key.column_iter().any(|col| col == *vec);

    // Check if the negative vector exists as a column in the matrix
    let as_column_neg = key.column_iter().any(|col| -col == *vec);

    as_column || as_column_neg
}

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
