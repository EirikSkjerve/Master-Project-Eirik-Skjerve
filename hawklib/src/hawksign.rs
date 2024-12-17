// all parameters
use crate::{
    hawkkeygen::gen_f_g,
    parameters::{hawk1024_params, hawk256_params, hawk512_params},
    rngcontext::{get_random_bytes, shake256x4, RngContext},
};

use crate::utils::*;

use sha3::{
    digest::{ExtendableOutputReset, Update},
    Shake256,
};

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

    res
}

// symbreak
pub fn symbreak(v: &Vec<i64>) -> bool {
    for x in v.iter() {
        if *x != 0 {
            if *x > 0 {
                return true;
            } else {
                return false;
            }
        }
    }
    false
}

pub fn sample(s: &[u8], t: Vec<u8>, n: usize) -> Vec<i64> {
    // get the correct tables
    let (t0, t1) = match n {
        256 => (hawk256_params::T0, hawk256_params::T1),
        512 => (hawk512_params::T0, hawk512_params::T1),
        _ => (hawk1024_params::T0, hawk1024_params::T1),
    };

    // vector y of random high numbers
    // note that the entries in y are uniformly distributed

    let y = shake256x4(s, 5 * n / 2);

    // the value that is used as a part of the sample
    let mut v: i64;
    // initialize empty vector for the sample vector
    let mut x: Vec<i64> = vec![0; 2 * n];

    // since y is the result of 4 interleaved shake256 instances
    // the following indexing will access them in an appropriate manner
    for j in 0..4 {
        for i in 0..(n / 8) {
            for k in 0..4 {
                // this is the current index for our point x
                let r = 16 * i + 4 * j + k;

                let a = y[(j + 4 * ((5 * i) + k)) as usize];

                let b = modulo(y[j + 4 * (5 * i + 4)] >> (16 * k), 1 << 15);

                // the final scaled number we are using
                let c = modulo(a as u128, 1 << 63) + (1 << 63) * b as u128;

                // initialize v0, v1, and z to zero
                let mut v0: i64 = 0;
                let mut v1: i64 = 0;
                let mut z = 0;

                // here the actual sampling is done by checking how many of the elements in the
                // CD-Tables are strictly larger than the uniformly sampled c
                // and since the tables are CDTs of discrete gaussian, the outputted number is a
                // random variable from the discrete gaussian distribution
                loop {
                    if t0[z] == 0 || t1[z] == 0 {
                        break;
                    }
                    if c < t0[z] {
                        v0 += 1;
                    }
                    if c < t1[z] {
                        v1 += 1;
                    }
                    z += 1;
                }

                // check the target vector at the current index
                // and add v0 or v1 based on this
                if t[r as usize] == 0 {
                    v = 2 * v0;
                } else {
                    v = 2 * v1 + 1;
                }

                // flip the sign if the original value from y is too high
                if a >= 1 << 63 {
                    // println!("{}", r);
                    v = -v;
                }

                // add the sample to vector x
                x[r as usize] = v;
            }
        }
    }
    // return x
    x
}

pub fn hawksign(
    privkey: &(Vec<u8>, Vec<i64>, Vec<i64>),
    msg: &[u8],
    n: usize,
) -> (Vec<i64>, Vec<u8>) {
    //
    // given secret key components and message, compute a signature
    //

    let (kgseed, bigf, bigg) = privkey;

    assert_eq!(bigf.len(), n);
    assert_eq!(bigg.len(), n);
    assert!(n == 256 || n == 512 || n == 1024);

    // get the right parameters
    let (lensalt, sigmaverify) = match n {
        256 => (hawk256_params::LENSALT, hawk256_params::SIGMAVERIFY),
        512 => (hawk512_params::LENSALT, hawk512_params::SIGMAVERIFY),
        _ => (hawk1024_params::LENSALT, hawk1024_params::SIGMAVERIFY),
    };

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

    // F and G mod 2
    let bigf_mod2: Vec<i64> = poly_mod2(&bigf);
    let bigg_mod2: Vec<i64> = poly_mod2(&bigg);

    // start a loop that terminates when a valid signature is created
    loop {
        a += 2;

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
        let t0 = poly_add(
            &poly_mult_ntt(&h0, &f, p),
            &poly_mult_ntt(&h1, &bigf_mod2, p),
        );

        let t1 = poly_add(
            &poly_mult_ntt(&h0, &g, p),
            &poly_mult_ntt(&h1, &bigg_mod2, p),
        );

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
        let x = sample(&s, t, n);

        // split x into two vectors
        let (x0, x1) = (&x[0..n].to_vec(), &x[n..].to_vec());

        // check norm of vector x is not too high
        // bounded by 8n*(sigma^2)
        if (l2norm(&x) as f64) > (8 * n) as f64 * sigmaverify.powi(2) {
            continue;
        }

        // compute one part of the signature
        // w = B^-1 x, so w1 = g*x0 - f*x1
        // let mut w0 = poly_sub(&poly_mult_ntt(&bigg, &x0, p), &poly_mult_ntt(&bigf, &x1, p));
        let mut w1 = poly_sub(&poly_mult_ntt(&f, &x1, p), &poly_mult_ntt(&g, &x0, p));

        // let mut w0_m = poly_sub(&polymul(&bigg, &x0), &polymul(&bigf, &x1));
        // let mut w1_m = poly_sub(&polymul(&f, &x1), &polymul(&g, &x0));
        //
        // assert_eq!(w0, w0_m);
        // assert_eq!(w1, w1_m);
        //
        // let rotf = rot(&f);
        // let rotg = rot(&g);
        // let rotbigf = rot(&bigf);
        // let rotbigg = rot(&bigg);
        //
        // let w0_2 = poly_sub(&matmul(&rotbigg, &x0), &matmul(&rotbigf, &x1));
        // let w1_2 = poly_sub(&matmul(&rotf, &x1), &matmul(&rotg, &x0));
        // assert_eq!(w0, w0_2);
        // assert_eq!(w1, w1_2);
        //
        // let rotbinv = rot_key(
        //     &bigg,
        //     &g.clone().iter().map(|&x| -x).collect(),
        //     &bigf.clone().iter().map(|&x| -x).collect(),
        //     &f,
        // );
        //
        // let w_t = matmul(&rotbinv, &x);
        //
        // let mut w_t_c = w0.clone();
        // w_t_c.append(&mut w1.clone());
        // assert_eq!(w_t_c, w_t);

        // check symbreak condition
        if !symbreak(&w1) {
            // w0 = w0.iter().map(|&x| -x).collect();
            w1 = w1.iter().map(|&x| -x).collect();
        }
        // println!("w0 from sign: {:?}", w0);

        // compute the actual signature sig = (h-w)/2
        let sig: Vec<i64> = poly_sub(&h1, &w1).iter().map(|&x| x >> 1).collect();
        // println!("s: {:?}", sig);

        return (sig, salt);
    }
}

pub fn hawksign_total(
    privkey: &(Vec<u8>, Vec<i64>, Vec<i64>),
    msg: &[u8],
    n: usize,
) -> (Vec<i64>, Vec<i64>) {
    //
    // given secret key components and message, compute a signature
    // unlike specifications, return the entire w as signature
    // so that one does not need to use "rebuild" function to get entire w
    // also returns the raw x-vector

    let (kgseed, bigf, bigg) = privkey;

    assert_eq!(bigf.len(), n);
    assert_eq!(bigg.len(), n);
    assert!(n == 256 || n == 512 || n == 1024);

    // get the right parameters
    let (lensalt, sigmaverify) = match n {
        256 => (hawk256_params::LENSALT, hawk256_params::SIGMAVERIFY),
        512 => (hawk512_params::LENSALT, hawk512_params::SIGMAVERIFY),
        _ => (hawk1024_params::LENSALT, hawk1024_params::SIGMAVERIFY),
    };

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

    // F and G mod 2
    let bigf_mod2: Vec<i64> = poly_mod2(&bigf);
    let bigg_mod2: Vec<i64> = poly_mod2(&bigg);

    // start a loop that terminates when a valid signature is created
    loop {
        a += 2;

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
        let t0 = poly_add(
            &poly_mult_ntt(&h0, &f, p),
            &poly_mult_ntt(&h1, &bigf_mod2, p),
        );

        let t1 = poly_add(
            &poly_mult_ntt(&h0, &g, p),
            &poly_mult_ntt(&h1, &bigg_mod2, p),
        );

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
            rng.random(80).to_vec(),
        ]);

        // sample vector x = (x0, x1)
        let x = sample(&s, t, n);

        // split x into two vectors
        let (x0, x1) = (&x[0..n].to_vec(), &x[n..].to_vec());

        // check norm of vector x is not too high
        // bounded by 8n*(sigma^2)
        if (l2norm(&x) as f64) > (8 * n) as f64 * sigmaverify.powi(2) {
            continue;
        }

        // compute one part of the signature
        // w = B^-1 x, so w1 = g*x0 - f*x1
        let mut w0 = poly_sub(&poly_mult_ntt(&bigg, &x0, p), &poly_mult_ntt(&bigf, &x1, p));
        let mut w1 = poly_sub(&poly_mult_ntt(&f, &x1, p), &poly_mult_ntt(&g, &x0, p));

        // check symbreak condition
        if !symbreak(&w1) {
            w0 = w0.iter().map(|&x| -x).collect();
            w1 = w1.iter().map(|&x| -x).collect();
        }

        // in this version, return only w, not s=(h-w)/2
        let mut w = w0.clone();
        w.append(&mut w1.clone());
        return (w, x);
    }
}


pub fn hawksign_x_only(
    privkey: &(Vec<u8>, Vec<i64>, Vec<i64>),
    msg: &[u8],
    n: usize,
) -> Vec<i64> {
    //
    // given secret key components and message, compute a signature
    // unlike specifications, return the entire w as signature
    // so that one does not need to use "rebuild" function to get entire w
    // also returns the raw x-vector

    let (kgseed, bigf, bigg) = privkey;

    assert_eq!(bigf.len(), n);
    assert_eq!(bigg.len(), n);
    assert!(n == 256 || n == 512 || n == 1024);

    // get the right parameters
    let (lensalt, sigmaverify) = match n {
        256 => (hawk256_params::LENSALT, hawk256_params::SIGMAVERIFY),
        512 => (hawk512_params::LENSALT, hawk512_params::SIGMAVERIFY),
        _ => (hawk1024_params::LENSALT, hawk1024_params::SIGMAVERIFY),
    };

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

    // F and G mod 2
    let bigf_mod2: Vec<i64> = poly_mod2(&bigf);
    let bigg_mod2: Vec<i64> = poly_mod2(&bigg);

    // start a loop that terminates when a valid signature is created
    loop {
        a += 2;

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
        let t0 = poly_add(
            &poly_mult_ntt(&h0, &f, p),
            &poly_mult_ntt(&h1, &bigf_mod2, p),
        );

        let t1 = poly_add(
            &poly_mult_ntt(&h0, &g, p),
            &poly_mult_ntt(&h1, &bigg_mod2, p),
        );

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
            rng.random(80).to_vec(),
        ]);

        // sample vector x = (x0, x1)
        let x = sample(&s, t, n);

        // check norm of vector x is not too high
        // bounded by 8n*(sigma^2)
        if (l2norm(&x) as f64) > (8 * n) as f64 * sigmaverify.powi(2) {
            continue;
        }

        return x;
    }
}
