// all parameters
use crate::{
    fft::inverse_fft,
    ntru_solve::ntrusolve,
    ntt_constants::res_z,
    parameters::{hawk1024_params, hawk256_params, hawk512_params},
};

use crate::utils::*;

// rng class
use crate::rngcontext::{get_random_bytes, shake256x4, RngContext};

//
// for simplicity I try and use only i64 for integers, BigInt for big integers, and f64 for floats
//

// inner computations for hawk keygen
fn hawkkeygen_inner(
    n: usize,
    lenkgseed: usize,
    sigmakrsec: f64,
    beta0: f64,
    high11: usize,
    rng: &mut RngContext,
) -> Option<((Vec<u8>, Vec<i64>, Vec<i64>), (Vec<i64>, Vec<i64>))> {
    // create a seed that determines f and g for later regeneration
    let kgseed = rng.random(lenkgseed);

    // generate polynomials f and g
    let (f, g) = gen_f_g(&kgseed, n);

    // check that f and g have odd parity
    // I think this is mostly for encoding/decoding stuff, but just do it anyways

    if n == 256 || n == 512 || n == 1024 {

        if !(is_invertible(&f, 2) && is_invertible(&g, 2)) {
            // println!("Is not invertible");
            return None;
        }

        // check l2norm of f and g
        if ((l2norm(&f) + l2norm(&g)) as f64) <= 2.0 * (n as f64) * sigmakrsec.powi(2) {
            // println!("Too low norm");
            return None;
        }
    }

    // compute f* and g*
    let fadj = adjoint(&f);
    let gadj = adjoint(&g);

    // first prime for NTT computations
    let p = (1 << 16) + 1;

    // compute q00, first element of public key Q
    let q00 = poly_add(&poly_mult_ntt(&f, &fadj, p), &poly_mult_ntt(&g, &gadj, p));

    // prime numbers used for ntt computations with Q later
    let (p1, p2): (i64, i64) = (2147473409, 2147389441);

    // check q00 invertibility mod p1 and p2
    if !(is_invertible(&q00, p1) && is_invertible(&q00, p2)) {
        // println!("not invertible...");
        return None;
    }

    // compute inverse of q00 using fft
    let invq00 = inverse_fft(&q00);

    // check the constant term of 1/q00 is not too high
    if invq00[0] >= beta0 {
        // println!("Too high thing here");
        return None;
    }

    // calculate F and G as solutions to NTRU-equation
    // if there is no solution, function should be called again

    // let ntru_res: Option<(Vec<i64>, Vec<i64>)> = Some(ntrusolve(&f, &g));
    let ntrusolve_res = ntrusolve(&f, &g);
    match ntrusolve_res {
        Some(_) => {}
        None => {println!("Ntrusolve failed"); return None},
    }

    let (bigf, bigg) = ntrusolve_res.unwrap();

    // bound largest element on F and G
    if infnorm(&bigf) > 127 || infnorm(&bigg) > 127 {
        return None;
    }

    // compute F* and G*
    let bigfadj = adjoint(&bigf);
    let biggadj = adjoint(&bigg);

    // compute two more elements in Q
    let q01 = poly_add(
        &poly_mult_ntt(&bigf, &fadj, p),
        &poly_mult_ntt(&bigg, &gadj, p),
    );

    let p = 8380417;
    let q11 = poly_add(
        &poly_mult_ntt(&bigf, &bigfadj, p),
        &poly_mult_ntt(&bigg, &biggadj, p),
    );

    // set bound on coefficients of q11
    for i in 1..q11.len() {
        if q11[i].abs() >= 1 << high11 {
            return None;
        }
    }

    // private key is seed for generating f and g, and F and G
    let private_key = (kgseed, bigf.clone(), bigg.clone());

    // pubblic key is two elements from Q. The rest can be reconstructed
    let public_key = (q00, q01);

    // note that we do not perform compression here

    // return generated keypairs
    Some((private_key, public_key))
}

pub fn hawkkeygen(n: usize, seed: Option<Vec<u8>>) -> ((Vec<u8>, Vec<i64>, Vec<i64>), (Vec<i64>, Vec<i64>)) {
    // main function for generating hawk keypair
    // input degree n
    // outputs private key (kgseed, F, G)
    //         public key  (q00, q01)

    let lenkgseed = match n {
        256 => hawk256_params::LENKGSEED,
        512 => hawk512_params::LENKGSEED,
        1024 => hawk1024_params::LENKGSEED,
        _ => hawk256_params::LENKGSEED,
    };

    let sigmakrsec = match n {
        256 => hawk256_params::SIGMAKRSEC,
        512 => hawk512_params::SIGMAKRSEC,
        1024 => hawk1024_params::SIGMAKRSEC,
        _ => hawk256_params::SIGMAKRSEC,
    };

    let beta0 = match n {
        256 => hawk256_params::BETA0,
        512 => hawk512_params::BETA0,
        1024 => hawk1024_params::BETA0,
        _ => hawk256_params::BETA0,
    };

    let high11 = match n {
        256 => hawk256_params::HIGH11,
        512 => hawk512_params::HIGH11,
        1024 => hawk1024_params::HIGH11,
        _ => hawk256_params::HIGH11,
    };

    // initialize a new RngContext instance, used for generating random bits
    let mut rng = RngContext::new(&seed.unwrap_or(get_random_bytes(100)));

    // reset precomputed roots for usage in ntt
    res_z();

    // try to generate hawk keypair
    // if inner function fails, loop runs again
    loop {
        match hawkkeygen_inner(n, lenkgseed, sigmakrsec, beta0, high11, &mut rng) {
            Some(keypair) => return keypair,
            None => {}
        }
    }
}

// generates polynomials f and g
pub fn gen_f_g(seed: &[u8], n: usize) -> (Vec<i64>, Vec<i64>) {
    let b = n / 64;

    // get an array of 64-bit values
    let y = shake256x4(&seed, 2 * n * b / 64);

    // construct a sequence of bits from y
    let mut ybits: Vec<u8> = vec![0; b * 2 * n];
    for (j, y) in y.iter().enumerate() {
        for bi in 0..64 {
            ybits[j * 64 + bi] = ((y >> bi) & 1) as u8;
        }
    }

    // generate f and g from centered binomial distribution
    // if e.g. n = 256, b = 4, so f and g consists of random numbes from the interval [-2,-1,0,1,2]
    // if n = 512, b = 8, interval = [-4,...,4]
    // if n = 1024, b = 16, interval = [-8,...,8]
    let mut f: Vec<i64> = vec![0; n];
    let mut sum;

    // get some bounded number from ybits
    for i in 0..n {
        sum = 0;
        for j in 0..b {
            sum += ybits[i * b + j] as i64;
        }

        // center the number around 0
        f[i] = sum - (b / 2) as i64;
    }

    // reuse "sum" variable here
    let mut g: Vec<i64> = vec![0; n];

    for i in 0..n {
        sum = 0;
        for j in 0..b {
            sum += ybits[(i + n) * b + j] as i64;
        }

        g[i] = sum - (b / 2) as i64;
    }

    // returns the two vectors
    return (f, g);
}
