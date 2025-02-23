use sha3::{
    digest::{ExtendableOutputReset, Update},
    Shake256,
};

use crate::hawksign::symbreak;
use crate::parameters::{hawk1024_params, hawk256_params, hawk512_params};
use crate::utils::{bytes_to_poly, modulo, poly_sub};
use crate::verifyutils::{poly_qnorm, rebuildw0};

// poly times const
pub fn ptc(a: &Vec<i64>, b: i64) -> Vec<i64> {
    a.clone().iter().map(|&x| x * b).collect()
}

pub fn hawkverify(
    msg: &[u8],
    pub_key: &(Vec<i64>, Vec<i64>),
    signature: &Vec<i64>,
    salt: &Vec<u8>,
    n: usize,
) -> bool {
    //
    // main method for verifying hawk signature on message
    //

    // get the correct parameters
    let (highs0, highs1, high00, high01, sigmaverify) = match n {
        256 => (
            hawk256_params::HIGHS0,
            hawk256_params::HIGHS1,
            hawk256_params::HIGH00,
            hawk256_params::HIGH01,
            hawk256_params::SIGMAVERIFY,
        ),

        512 => (
            hawk512_params::HIGHS0,
            hawk512_params::HIGHS1,
            hawk512_params::HIGH00,
            hawk512_params::HIGH01,
            hawk512_params::SIGMAVERIFY,
        ),

        1024 => (
            hawk1024_params::HIGHS0,
            hawk1024_params::HIGHS1,
            hawk1024_params::HIGH00,
            hawk1024_params::HIGH01,
            hawk1024_params::SIGMAVERIFY,
        ),

        _ => (
            hawk256_params::HIGHS0,
            hawk256_params::HIGHS1,
            hawk256_params::HIGH00,
            hawk256_params::HIGH01,
            hawk256_params::SIGMAVERIFY,
        ),
    };

    // convert signature to Vec<i64>
    let s1: Vec<i64> = signature.iter().map(|&x| x as i64).collect();
    let (q00, q01) = pub_key;

    // compute hash digest of message m
    let mut shaker = Shake256::default();
    shaker.update(msg);
    let mut m: Vec<u8> = vec![0; 64];
    shaker.finalize_xof_reset_into(&mut m);

    // create buffer for vector h whose length depends on degree n
    let mut h: Vec<u8> = vec![0; n / 4];

    // digest h is digest of message+salt
    shaker.update(&m);
    shaker.update(&salt);

    // compute digest
    shaker.finalize_xof_reset_into(&mut h);

    // convert digest h to usable polynomials
    let (h0, h1) = (
        &bytes_to_poly(&h[0..n / 8], n),
        &bytes_to_poly(&h[n / 8..n / 4], n),
    );

    // reconstruct digest h

    let w1 = poly_sub(&h1, &ptc(&s1, 2));

    if !symbreak(&w1) {
        println!("Symbreak failed");
        return false;
    }

    // rebuild the missing part of the signature w
    let w0 = rebuildw0(&q00, &q01, &w1, &h0, highs0, highs1, high00, high01);

    // primes used for doing ntt computations with Q
    let (p1, p2): (i64, i64) = (2147473409, 2147389441);

    // calculate q-norm of signature w.r.t. two primes p1 and p2, same as used in signature
    // generation
    let r1 = poly_qnorm(&q00, &q01, &w0, &w1, p1);
    let r2 = poly_qnorm(&q00, &q01, &w0, &w1, p2);

    if r1 != r2 || modulo(r1, n as i64) != 0 {
        println!("r1 != r2. Returning false");
        return false;
    }

    let r1 = r1 / (n as i64);

    // check the q norm is not too high
    if (r1 as f64) > (8 * n) as f64 * sigmaverify.powi(2) {
        return false;
    }

    // accept signature
    true
}
