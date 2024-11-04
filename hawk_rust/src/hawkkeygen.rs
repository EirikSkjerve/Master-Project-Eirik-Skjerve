// all parameters
use crate::parameters::{hawk1024_params, hawk256_params, hawk512_params};

// shake256
use sha3::{
    digest::{ExtendableOutput, Update},
    Shake256,
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
) {
    // create a seed that determines f and g for later regeneration
    let kgseed = rng.random(lenkgseed);

    // generate polynomials f and g
    let (f, g) = gen_f_g(&kgseed, n);

    // check that f and g have odd parity
    // I think this is mostly for encoding/decoding stuff, but just do it anyways
    if !(is_invertible(&f, 2) && is_invertible(&g, 2)) {
        return;
    }

    // check l2norm of f and g
    if ((l2norm(&f) + l2norm(&g)) as f64) <= 2.0 * (n as f64) * sigmakrsec.powi(2) {
        return;
    }

    // compute f* and g*
    let fadj = adjoint(&f);
    let gadj = adjoint(&g);

    // first prime for NTT computations
    let p = (1 << 16) + 1;

    // compute q00, first element of public key Q
    let q00 = poly_add(&poly_mult_ntt(&f, &fadj, p), &poly_mult_ntt(&g, &gadj, p));
}

pub fn hawkkeygen_256() {
    const N: usize = 256;

    // initialize a new RngContext instance, used for generating random bits
    let mut rng = RngContext::new(&get_random_bytes(10));

    loop {
        hawkkeygen_inner(
            N,
            hawk256_params::LENKGSEED,
            hawk256_params::SIGMAKRSEC,
            hawk256_params::BETA0,
            hawk256_params::HIGH11,
            &mut rng,
        );
    }
}

pub fn hawkkeygen_512() {
    const N: usize = 512;
}

pub fn hawkkeygen_1024() {
    const N: usize = 1024;
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
