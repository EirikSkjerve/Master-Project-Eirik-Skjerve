use crate::fft::{inverse_fft, mul_fft_i64};
use crate::ntru_solve::ntrusolve;
use crate::rngcontext::{shake256x4, RngContext};
use crate::utils::{
    adjoint, bigint_to_i64_vec, bigint_vec, infnorm, is_invertible, l2norm, poly_add, poly_mult_ntt,
};
use crate::codec::enc_pub;

use num_bigint::{BigInt, BigUint, ToBigInt, ToBigUint};
use num_traits::{One, Signed, ToPrimitive, Zero};

pub fn hawkkeygen(
    logn: u8,
    initial_seed: usize,
) -> (
    Vec<i64>,
    Vec<i64>,
    Vec<i64>,
    Vec<i64>,
    Vec<i64>,
    Vec<i64>,
    usize,
    i32,
) {
    // not doing recursion
    let mut rng = RngContext::new(initial_seed as u128);
    let mut counter = 0;
    loop {
        counter += 1;
        // for each new loop, kgseed will be a new random value
        let kgseed = rng.rnd(128) as usize;

        // generate f and g from centered binomial distribution
        let f_g = generate_f_g(kgseed, logn);
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

        // should have some test if ntrusolve fails
        let (F, G) = ntrusolve(bigint_vec(&f), bigint_vec(&g));

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
        let q01 = poly_add(
            &poly_mult_ntt(&F, &F_star, p),
            &poly_mult_ntt(&G, &G_star, p),
        );

        let encoded = enc_pub(logn as usize,&q00, &q01);
        println!("encoded: {:?}", encoded);

        return (f.clone(), g.clone(), F, G, q00, q01, kgseed, counter);
    }
}

// generates polynomials f and g
pub fn generate_f_g(seed: usize, logn: u8) -> (Vec<i64>, Vec<i64>) {
    // expand logn -> n
    let n = 1 << logn;
    let b = n / 64;
    assert!(b == 4 || b == 8 || b == 16);

    // get an array of 64-bit values
    let y = shake256x4(&seed.to_ne_bytes(), 2 * n * b / 64);

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
