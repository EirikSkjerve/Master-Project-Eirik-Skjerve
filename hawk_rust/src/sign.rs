use sha3::{
    digest::{ExtendableOutput, ExtendableOutputReset, Update},
    Shake256,
};

use rand::prelude::*;
use rand::{rngs::StdRng, Rng, SeedableRng};

use crate::cdt::get_table;
use crate::keygen::generate_f_g;
use crate::rngcontext::{shake256x4, RngContext};
use crate::utils::{bytes_to_poly, modulo, poly_add, poly_mult_ntt, poly_sub};

use crate::codec::{dec_priv, enc_sig};

pub fn sample(seed: &[u8], t: Vec<u8>, n: usize) -> Vec<i8> {
    let (T0, T1) = get_table();
    let y = shake256x4(seed, 5 * n / 2);

    // following HAWK's implementation

    let mut v = 0;
    let mut x: Vec<i8> = vec![0; 2 * n];

    let base: u64 = 2;
    let base_pow_15 = base.pow(15);
    let base_pow_63 = base.pow(63) as u128;
    for j in 0..4 {
        for i in 0..(n / 8) {
            for k in 0..4 {
                let r = 16 * i + 4 * j + k;
                let a = y[(j + 4 * ((5 * i) + k)) as usize];
                let pow_temp = base.pow(16 * k as u32);

                let temp1 = y[(j + 4 * ((5 * i) + 4)) as usize] / pow_temp;
                let b = modulo(temp1, base_pow_15) as u128;

                let c = modulo(a as u128, base_pow_63) + base_pow_63 * b;

                let mut v0: i8 = 0;
                let mut v1: i8 = 0;
                let mut z = 0;

                loop {
                    if T0[z] == 0 || T1[z] == 0 {
                        break;
                    }
                    if c < T0[z] {
                        v0 += 1;
                    }
                    if c < T1[z] {
                        v1 += 1;
                    }
                    z += 1;
                }

                if t[r as usize] == 0 {
                    v = 2 * v0;
                } else {
                    v = 2 * v1 + 1;
                }

                if a >= base.pow(63) {
                    v = -v;
                }

                x[r as usize] = v;
            }
        }
    }
    return x;
}

pub fn sign(logn: usize, pk: &Vec<u8>, msg: usize) -> Vec<u8> {

    let (kgseed, Fmod2, Gmod2, hpub) = dec_priv(logn, pk);
    // this should be from logn
    // initialize a new RngContext with some random seed
    // use random() instead of fixed seed
    let mut seed_rng = rand::thread_rng();
    let mut rng = RngContext::new(seed_rng.gen());

    let n = 1 << logn;
    let (f, g) = generate_f_g(kgseed, logn);

    // compute hash M
    let mut shaker = Shake256::default();
    shaker.update(&msg.to_ne_bytes());
    let mut m: [u8; 64] = [0; 64];
    shaker.finalize_xof_reset_into(&mut m);

    let mut a: usize = 0;
    let p = (1 << 16) + 1;

    loop {
        // compute salt
        shaker.update(&m);
        shaker.update(&kgseed.to_ne_bytes());
        shaker.update(&a.to_ne_bytes());
        shaker.update(&rng.rnd(14).to_ne_bytes());
        // saltlen should be from parameters
        let mut salt: [u8; 14] = [0; 14];
        // resets the hasher instance
        shaker.finalize_xof_reset_into(&mut salt);

        // compute new hash h

        shaker.update(&m);
        shaker.update(&salt);
        // the 256 is the degree. Should depend on input logn
        let mut h: [u8; 256 / 4] = [0; 256 / 4];
        shaker.finalize_xof_reset_into(&mut h);

        // convert h to two polynomials
        let (h0, h1) = (
            &bytes_to_poly(&h[0..256 / 8], n),
            &bytes_to_poly(&h[(256 / 8)..256 / 4], n),
        );



        // compute target vectors t0, t1

        let mut t0: Vec<u8> = Vec::with_capacity(n);
        let mut t1: Vec<u8> = Vec::with_capacity(n);

        let temp_t0 = poly_add(&poly_mult_ntt(&h0, &f, p), &poly_mult_ntt(&h1, &Fmod2, p));
        let temp_t1 = poly_add(&poly_mult_ntt(&h0, &g, p), &poly_mult_ntt(&h1, &Gmod2, p));

        for i in 0..n {
            // we can be sure these values fit inside an u8 since values are 0 and 1
            t0.push(modulo(temp_t0[i], 2) as u8);
            t1.push(modulo(temp_t1[i], 2) as u8);
        }

        let t = concat_bytes(&vec![t0, t1]);
        // get random seed = M || kgseed || a+1 || rnd(320)
        let seed = rng.rnd(40);
        let kgseed_b = to_bytes_sized(kgseed, 64);
        let arr = vec![
            m.to_vec(),
            kgseed_b,
            (a + 1).to_ne_bytes().to_vec(),
            seed.to_ne_bytes().to_vec(),
        ];
        // this should be concatenation
        let s_temp = concat_bytes(&arr);
        let s = vec_to_slice(&s_temp);

        // compute (x0, x1) from sample()

        let x = sample(s, t, n);

        let x0 = &x[0..n];
        let x1 = &x[n..];

        // increment for new salt if failure
        a += 2;

        let sigmaverify: f64 = 1.024;
        let factor: f64 = (8 * n) as f64;
        let l2normsum = l2norm_sign(x0) + l2norm_sign(x1);
        println!("l2norm(x0) + l2norm(x1) = {}", l2normsum);
        // continue loop if some requirements are not fulfilled
        if (l2normsum as f64) > factor * sigmaverify.powi(2) {
            continue;
        }

        // convert x0 and x1 to Vec<i64>

        let x0_i64: Vec<i64> = x0.iter().map(|&x| x as i64).collect();
        let x1_i64: Vec<i64> = x1.iter().map(|&x| x as i64).collect();

        let mut w1 = poly_sub(
            &poly_mult_ntt(&f, &x1_i64, p),
            &poly_mult_ntt(&g, &x0_i64, p),
        );

        if !symbreak(&w1) {
            w1 = w1.iter().map(|&x| -x).collect();
        }

        println!("w1 from sign: {:?}", w1);
        let sig: Vec<i64> = poly_sub(&h1, &w1).iter().map(|&x| x>>1).collect();

        // println!("sig un-encoded from sig: {:?}", sig);

        let sig_enc = enc_sig(logn, &salt.to_vec(), &sig);
        if sig_enc[0] == 0 && sig_enc.len() == 1 {
            continue
        }
        println!("s1 = {:?}", sig);
        println!("h0 = {:?} \nh1 = {:?}", h0, h1);
        return sig_enc;
    }
}

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
    // converts to f64 for usage in signing procedure
    let mut sum: i64 = 0;
    for i in 0..a.len() {
        sum += (a[i]*a[i]) as i64;
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

fn to_bytes_sized(a: usize, size: usize) -> Vec<u8> {
    let a_b = a.to_ne_bytes();
    assert!(a_b.len() < size);
    let mut res = Vec::with_capacity(size);
    for i in 0..size {
        if i < a_b.len() {
            res.push(a_b[i]);
        } else {
            res.push(0);
        }
    }
    return res;
}

fn add_bytes(arr: Vec<Vec<u8>>) -> Vec<u64> {
    /*
     * adds together an arbitrary number of byte arrays into a vector
     */
    // find the max size array
    let mut n_temp = 0;
    for a in arr.iter() {
        let a_max = a.iter().max().unwrap();
        if *a_max > n_temp {
            n_temp = *a_max;
        }
    }

    let n = n_temp as usize;

    // initialize vector
    let mut res: Vec<u64> = Vec::with_capacity(n);
    let mut temp = 0;
    for i in 0..n {
        // sum together all values at index i
        temp = 0;
        for a in arr.iter() {
            if i < a.len() {
                temp += a[i];
            }
        }
        res.push(temp as u64);
    }
    return res;
}
