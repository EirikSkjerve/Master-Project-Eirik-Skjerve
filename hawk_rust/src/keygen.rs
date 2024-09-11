use crate::codec::{enc_priv, enc_pub};
use crate::fft::{inverse_fft, mul_fft_i64};
use crate::ntru_solve::ntrusolve;
use crate::params::params_i;
use crate::rngcontext::{shake256x4, RngContext};
use crate::utils::{
    adjoint, bigint_to_i64_vec, bigint_vec, infnorm, is_invertible, l2norm, poly_add, poly_mult_ntt,
};

use sha3::{
    digest::{ExtendableOutput, ExtendableOutputReset, Update},
    Shake256,
};

use num_bigint::{BigInt, BigUint, ToBigInt, ToBigUint};
use num_traits::{One, Signed, ToPrimitive, Zero};

pub fn hawkkeygen(logn: usize, initial_seed: &[u8]) -> (Vec<u8>, Vec<u8>) {
    // not doing recursion
    let mut rng = RngContext::new(initial_seed);
    let mut counter = 0;
    loop {
        counter += 1;
        // for each new loop, kgseed will be a new random value
        let kgseed = rng.random(128 / 8);

        // generate f and g from centered binomial distribution
        let f_g = generate_f_g(&kgseed, logn);
        let f = f_g.0;
        let g = f_g.1;

        if !is_invertible(&f, 2) || !is_invertible(&g, 2) {
            // println!("Restart 1");
            continue;
        }

        let n = 1 << logn;

        if ((l2norm(&f) + l2norm(&g)) as f64) <= (2.0 * (n as f64) * (1.042 as f64).powi(2)) {
            // println!("Restart 2");
            continue;
        }

        let f_star = adjoint(&f);
        let g_star = adjoint(&g);

        let p = (1 << 16) + 1;

        let q00 = poly_add(
            &poly_mult_ntt(&f, &f_star, p),
            &poly_mult_ntt(&g, &g_star, p),
        );

        let p1: u32 = 2147473409;
        let p2: u32 = 2147389441;

        if !is_invertible(&q00, p1) || !is_invertible(&q00, p2) {
            // println!("Restart 3");
            continue;
        }

        let invq00 = inverse_fft(&q00);

        if invq00[0] >= 1.0 / 250.0 {
            // println!("Restart 4");
            continue;
        }

        // calculate F and G
        let (F, G) = ntrusolve(bigint_vec(&f), bigint_vec(&g));

        // if F and G are not found, retry
        if (F.len() == 1 && F[0] == BigInt::ZERO) && (G.len() == 1 && G[0] == BigInt::ZERO) {
            println!("Could not calculate F and G. Retrying");
            continue;
        }

        let (F, G) = (bigint_to_i64_vec(F), bigint_to_i64_vec(G));

        if infnorm(&F) > 127 || infnorm(&G) > 127 {
            println!("Restart 5");
            continue;
        }

        let F_star = adjoint(&F);
        let G_star = adjoint(&G);

        let q01 = poly_add(
            &poly_mult_ntt(&F, &f_star, p),
            &poly_mult_ntt(&G, &g_star, p),
        );

        let p = 8380417;
        let q11 = poly_add(
            &poly_mult_ntt(&F, &F_star, p),
            &poly_mult_ntt(&G, &G_star, p),
        );

        let mut flag = false;
        for i in 1..q11.len() {
            if q11[i].abs() >= 1 << (params_i(logn as usize, "high11")) {
                flag = true;
            }
        }
        if flag {
            continue;
        }

        let pub_enc = enc_pub(logn as usize, &q00, &q01);

        // encoding failed if the following is true
        if pub_enc[0] == 0 && pub_enc.len() == 1 {
            continue;
        }

        let mut shaker = Shake256::default();
        let pk: &[u8] = &pub_enc;
        shaker.update(pk);
        // this should be retrieved from params later
        let mut hpub: [u8; 16] = [0; 16];
        shaker.finalize_xof_into(&mut hpub);

        let priv_enc = enc_priv(&kgseed, &F, &G, &hpub);
        // encode private key as encode_private(kgseed, F, G, hpub);

        return (priv_enc, pub_enc);
    }
}

// generates polynomials f and g
pub fn generate_f_g(seed: &[u8], logn: usize) -> (Vec<i64>, Vec<i64>) {
    // expand logn -> n
    let n = 1 << logn;
    let b = n / 64;
    assert!(b == 4 || b == 8 || b == 16);

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
