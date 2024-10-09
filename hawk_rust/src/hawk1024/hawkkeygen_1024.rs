// all utils files
use crate::hawk1024::codec_1024::{enc_priv, enc_pub};
use crate::fft::inverse_fft;
use crate::ntru_solve::ntrusolve;
use crate::rngcontext::{shake256x4, RngContext};
use crate::utils::{
    adjoint, bigint_to_i64_vec, bigint_vec, infnorm, is_invertible, l2norm, poly_add, poly_mult_ntt,
};

// shake256 
use sha3::{
    digest::{ExtendableOutput, Update},
    Shake256,
};

// parameters for hawk 1024
use crate::parameters::hawk1024_params::*;
use num_bigint::BigInt;


/// Generates a HAWK 512 public/private key pair
/// Will return encoded/compressed keys
pub fn hawkkeygen_1024(initial_seed: &[u8]) -> (Vec<u8>, Vec<u8>) {

    let logn = 10;
    // initialize a new RngContext instance, used for generating random bits
    let mut rng = RngContext::new(initial_seed);

    // start main loop 
    // loop might need to run several times
    loop {
        // for each new loop, kgseed will be a new random value
        let kgseed = rng.random(LENKGSEED);

        // generate f and g from a centered binomial distribution
        let f_g = generate_f_g(&kgseed, logn);
        let f = f_g.0;
        let g = f_g.1;

        // check invertibility mod 2
        if !is_invertible(&f, 2) || !is_invertible(&g, 2) {
            continue;
        }

        let n = 1 << logn;

        // check that norm of f and g is not too low
        if ((l2norm(&f) + l2norm(&g)) as f64) <= 2.0*(n as f64) * SIGMAKRSEC.powi(2) {
            continue;
        }

        // compute the Hermitian Adjoint of f and g
        let f_star = adjoint(&f);
        let g_star = adjoint(&g);

        // prime number used for ntt computations
        let p = (1 << 16) + 1;

        // compute first element of public key Q
        let q00 = poly_add(
            &poly_mult_ntt(&f, &f_star, p),
            &poly_mult_ntt(&g, &g_star, p),
        );

        // two prime numbers used for computation with Q as quadratic form
        let p1: u32 = 2147473409;
        let p2: u32 = 2147389441;

        // check invertibility of q00 mod these two primes
        if !is_invertible(&q00, p1) || !is_invertible(&q00, p2) {
            continue;
        }

        // compute the inverse of q00 using fft
        let invq00 = inverse_fft(&q00);

        // check the constant term of 1/q00 is not too high
        if invq00[0] >= BETA0 {
            continue;
        }

        // calculate bigf and bigg
        // here we convert f and g into vectors of the BigInt type
        let (bigf, bigg) = ntrusolve(bigint_vec(&f), bigint_vec(&g));

        // if bigf and bigg are not found, retry
        if (bigf.len() == 1 && bigf[0] == BigInt::ZERO) && (bigg.len() == 1 && bigg[0] == BigInt::ZERO) {
            continue;
        }

        // convert the output from ntrusolve (which are vectors of BigInt) into i64 form
        let (bigf, bigg) = (bigint_to_i64_vec(bigf), bigint_to_i64_vec(bigg));

        // bound on the largest element of bigf and bigg
        if infnorm(&bigf) > 127 || infnorm(&bigg) > 127 {
            continue;
        }

        // compute Hermitian Adjoint of F and G
        let bigf_star = adjoint(&bigf);
        let bigg_star = adjoint(&bigg);

        // compute two more elements in Q
        let q01 = poly_add(
            &poly_mult_ntt(&bigf, &f_star, p),
            &poly_mult_ntt(&bigg, &g_star, p),
        );

        let p = 8380417;
        let q11 = poly_add(
            &poly_mult_ntt(&bigf, &bigf_star, p),
            &poly_mult_ntt(&bigg, &bigg_star, p),
        );

        // set bound on coefficients of q11
        let mut flag = false;
        for i in 1..q11.len() {
            if q11[i].abs() >= 1 << HIGH11 {
                flag = true;
            }
        }
        if flag {
            continue;
        }

        // encode the public key
        let pub_enc = enc_pub(logn as usize, &q00, &q01);

        // encoding failed if the following is true
        if pub_enc[0] == 0 && pub_enc.len() == 1 {
            continue;
        }

        // compute a hash digest of the (encoded) public key
        let mut shaker = Shake256::default();

        let pk: &[u8] = &pub_enc;
        shaker.update(pk);

        let mut hpub: [u8; LENHPUB] = [0; LENHPUB];
        shaker.finalize_xof_into(&mut hpub);

        // encode private key
        let priv_enc = enc_priv(&kgseed, &bigf, &bigg, &hpub);

        // return the keys
        return (priv_enc, pub_enc);
    }
}

// generates polynomials f and g
pub fn generate_f_g(seed: &[u8], logn: usize) -> (Vec<i64>, Vec<i64>) {
    // expand logn -> n
    let n = 1 << logn;
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
