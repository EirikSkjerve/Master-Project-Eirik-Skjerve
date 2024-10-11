use sha3::{
    digest::{ExtendableOutputReset, Update},
    Shake256,
};

use crate::hawk1024::hawkkeygen_1024::generate_f_g;
use crate::rngcontext::{get_random_bytes, shake256x4, RngContext};
use crate::utils::{bytes_to_poly, modulo, poly_add, poly_mult_ntt, poly_sub};

use crate::hawk1024::codec_1024::{dec_priv, enc_sig};

use crate::parameters::hawk1024_params::*;

const N: usize = 1024;

pub fn sample(seed: &[u8], t: Vec<u8>) -> [i8; 2*N] {
    // get the CDT for this degree
    let (t0, t1) = (T0, T1);

    // vector y of high numbers
    // note that the entries in y are uniformly distributed
    let y = shake256x4(seed, 5 * N / 2);

    // the value that is used as a part of the sample
    let mut v = 0;
    // initialize empty vector for the sample
    // let mut x: Vec<i8> = vec![0; 2 * n];
    let mut x: [i8; 2*N] = [0; 2*N];

    // since y is the result of 4 interleaved shake256 instances
    // the following indexing will access them in an appropriate manner
    for j in 0..4 {
        for i in 0..(N / 8) {
            for k in 0..4 {
                // this is the current index for our point x
                let r = 16 * i + 4 * j + k;

                let a = y[(j + 4 * ((5 * i) + k)) as usize];

                let b = modulo(y[j + 4 * (5 * i + 4)] >> (16 * k), 1 << 15);

                // the final scaled number we are using
                let c = modulo(a as u128, 1 << 63) + (1 << 63) * b as u128;

                // initialize v0, v1, and z to zero
                let mut v0: i8 = 0;
                let mut v1: i8 = 0;
                let mut z = 0;

                // here the actual sampling is done by checking how many of the elements in the
                // CD-Tables are strictly larger than the uniformly sampled c
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
                // (higher than half???)
                if a >= 1 << 63 {
                    v = -v;
                }

                // add the sample to array x
                x[r as usize] = v;
            }
        }
    }
    return x;
}

pub fn hawksign_1024(sk: &Vec<u8>, msg: &[u8]) -> Vec<u8> {
    let logn = 10;
    let (kgseed, bigfmod2, biggmod2, hpub) = dec_priv(logn, sk);

    // convert the Vec<u8> kgseed to a &[u8]
    let kgseed = &kgseed;

    // initialize a new RngContext with some random seed
    // use random() instead of fixed seed
    // let mut rng = RngContext::new(seed_rng.gen());
    let mut rng = RngContext::new(&get_random_bytes(10));

    let (f, g) = generate_f_g(kgseed, logn);

    // compute hash M
    let mut shaker = Shake256::default();
    shaker.update(msg);
    shaker.update(&hpub[..]);
    let mut m: [u8; 64] = [0; 64];
    shaker.finalize_xof_reset_into(&mut m);

    let mut a: usize = 0;
    let p = (1 << 16) + 1;

    loop {
        // compute salt
        shaker.update(&m);
        shaker.update(&kgseed);
        shaker.update(&a.to_ne_bytes());
        shaker.update(&rng.random(LENSALT));
        // saltlen should be from parameters
        let mut salt: [u8; LENSALT] = [0; LENSALT];
        // resets the hasher instance
        shaker.finalize_xof_reset_into(&mut salt);

        // compute new hash h
        shaker.update(&m);
        shaker.update(&salt);

        let mut h: [u8; N / 4] = [0; N / 4];
        shaker.finalize_xof_reset_into(&mut h);

        // convert h to two polynomials
        let (h0, h1) = (
            &bytes_to_poly(&h[0..N / 8], N),
            &bytes_to_poly(&h[(N / 8)..N / 4], N),
        );

        // compute target vectors t0, t1
        // t = Bh
        let mut t0: Vec<u8> = Vec::with_capacity(N);
        let mut t1: Vec<u8> = Vec::with_capacity(N);

        let temp_t0 = poly_add(
            &poly_mult_ntt(&h0, &f, p),
            &poly_mult_ntt(&h1, &bigfmod2, p),
        );
        let temp_t1 = poly_add(
            &poly_mult_ntt(&h0, &g, p),
            &poly_mult_ntt(&h1, &biggmod2, p),
        );

        for i in 0..N {
            // we can be sure these values fit inside an u8 since values are 0 and 1
            t0.push(modulo(temp_t0[i], 2) as u8);
            t1.push(modulo(temp_t1[i], 2) as u8);
        }

        let t = concat_bytes(&vec![t0, t1]);

        // get random seed = M || kgseed || a+1 || rnd(320)
        let seed = rng.random(40);

        let arr = vec![
            m.to_vec(),
            kgseed.to_vec(),
            (a + 1).to_ne_bytes().to_vec(),
            seed.to_vec(),
        ];

        let s_temp = concat_bytes(&arr);
        let s = vec_to_slice(&s_temp);

        // compute (x0, x1) from sample()

        let x = sample(s, t.clone());

        let x0 = &x[0..N];
        let x1 = &x[N..];

        // increment for new salt if failure
        a += 2;

        let factor: f64 = (8 * N) as f64;
        let l2normsum = l2norm_sign(x0) + l2norm_sign(x1);

        // continue loop if some requirements are not fulfilled
        if (l2normsum as f64) > factor * SIGMAVERIFY.powi(2) {
            continue;
        }

        // convert x0 and x1 to Vec<i64>

        let x0_i64: Vec<i64> = x0.iter().map(|&x| x as i64).collect();
        let x1_i64: Vec<i64> = x1.iter().map(|&x| x as i64).collect();

        // compute one part of the signature
        // the remaining part is computed in signature verification
        let mut w1 = poly_sub(
            &poly_mult_ntt(&f, &x1_i64, p),
            &poly_mult_ntt(&g, &x0_i64, p),
        );

        if !symbreak(&w1) {
            w1 = w1.iter().map(|&x| -x).collect();
        }

        // this is the signature
        let sig: Vec<i64> = poly_sub(&h1, &w1).iter().map(|&x| x >> 1).collect();

        // encode the signature
        let sig_enc = enc_sig(logn, &salt.to_vec(), &sig);

        // restart if encoding fails
        if sig_enc[0] == 0 && sig_enc.len() == 1 {
            continue;
        }

        // return encoded signature
        return sig_enc;
    }
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
    return false;
}

pub fn l2norm_sign(a: &[i8]) -> i64 {
    // returns the l2 norm of polynomial/vector f as f[0]^2 + f[1]^2 +..+ f[n]^2
    // converts to i64 for usage in signing procedure
    let mut sum: i64 = 0;
    for i in 0..a.len() {
        sum += (a[i] as i64) * (a[i] as i64);
    }
    return sum;
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

fn vec_to_slice(v: &Vec<u8>) -> &[u8] {
    return v;
}
