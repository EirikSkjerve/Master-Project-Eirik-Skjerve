// all parameters
use crate::{
    fft::inverse_fft,
    hawkkeygen::gen_f_g,
    ntru_solve::ntrusolve,
    parameters::{hawk1024_params, hawk256_params, hawk512_params},
    rngcontext::{get_random_bytes, shake256x4, RngContext},
};

use crate::utils::*;

use sha3::{
    digest::{ExtendableOutputReset, Update},
    Shake256,
};

// might need some smart solution for sampling for different degrees because of tables T0 T1
//
// fn sample_inner(seed: &[u8], t: Vec<u8>, n: usize) -> Vec<i8> {}
//
//

fn poly_mod2(a: &Vec<i64>) -> Vec<i64> {
    // returns a copy of a polynomial with each coefficient reduced mod 2
    a.clone().iter().map(|&x| modulo(x, 2)).collect()
}

fn concat_bytes(arr: &Vec<Vec<u8>>) -> Vec<u8> {
    /*
     * concatenates a vector of vectors into a single vector
     */

    let sum: usize = arr.iter().map(|slice| slice.len()).sum();

    let mut res = Vec::with_capacity(sum);

    for a in arr.iter() {
        for i in 0..a.len() {
            res.push(a[i]);
        }
    }

    return res;
}

fn hawksign_inner(
    kgseed: Vec<u8>,
    bigf: Vec<i64>,
    bigg: Vec<i64>,
    msg: &[u8],
    n: usize,
    lensalt: usize,
) -> Option<Vec<u8>> {
    // create new rng
    let mut rng = RngContext::new(&get_random_bytes(10));

    // regenerate part of secret key
    let (f, g) = gen_f_g(&kgseed, n);

    // compute hash of message
    let mut shaker = Shake256::default();
    shaker.update(msg);

    // fixed size hash digest buffer
    let mut m: [u8; 64] = [0; 64];
    // compute digest in m and reset the shaker instance
    shaker.finalize_xof_reset_into(&mut m);

    // counter variable that helps make unique salts
    let mut a: usize = 0;

    // prime used for ntt computations
    let p = (1 << 16) + 1;

    // precompute f, g, F and G mod 2 to save computations
    let f: Vec<i64> = poly_mod2(&f);
    let g: Vec<i64> = poly_mod2(&g);
    let bigf: Vec<i64> = poly_mod2(&bigf);
    let bigg: Vec<i64> = poly_mod2(&bigg);

    loop {
        // compute salt
        shaker.update(&m);
        shaker.update(&kgseed);
        shaker.update(&a.to_ne_bytes());
        shaker.update(&rng.random(14));

        // create salt buffer with length from parameter
        let mut salt: Vec<u8> = vec![0; lensalt];
        // create salt digest
        shaker.finalize_xof_reset_into(&mut salt);

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

        // compute target vector t as B*h mod 2
        let t0 = poly_add(&poly_mult_ntt(&h0, &f, p), &poly_mult_ntt(&h1, &bigf, p));
        let t1 = poly_add(&poly_mult_ntt(&h0, &g, p), &poly_mult_ntt(&h1, &bigg, p));

        // join t0 and t1 together as Vec<u8>
        let t = concat_bytes(&vec![
            poly_mod2(&t0).iter().map(|&x| x as u8).collect(),
            poly_mod2(&t1).iter().map(|&x| x as u8).collect(),
        ]);

        // create seed for sampling of vector x
        // M || kgseed || a+1 || rnd(320)
        // here 320 bits <=> 80 bytes
        let s = concat_bytes(&vec![
            m.to_vec(),
            kgseed.to_vec(),
            (a + 1).to_ne_bytes().to_vec(),
            rng.random(40).to_vec(),
        ]);

        // sample vector x = (x0, x1)

        break;
    }

    None
}

pub fn hawksign_256(kgseed: Vec<u8>, bigf: Vec<i64>, bigg: Vec<i64>, msg: &[u8]) {
    const N: usize = 256;
    assert_eq!(bigf.len(), N);
    assert_eq!(bigg.len(), N);

    // before calling hawksign_inner with parameters, some computations needs to be done
    // here so that the compiler knows how much space to allocate
    // or maybe not??

    hawksign_inner(kgseed, bigf, bigg, msg, N, hawk256_params::LENSALT);
}
