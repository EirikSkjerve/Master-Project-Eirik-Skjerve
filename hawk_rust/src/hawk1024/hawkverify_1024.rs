
use sha3::{
    digest::{ExtendableOutputReset, Update},
    Shake256,
};

use crate::hawk1024::codec_1024::{dec_pub, dec_sig};
use crate::hawk1024::verifyutils_1024::*;
use crate::utils::{bytes_to_poly, modulo, poly_sub};
use crate::parameters::hawk1024_params::*;

pub fn hawkverify_1024(msg: &[u8], pub_key: &Vec<u8>, signature: &Vec<u8>) -> bool {
    let logn = 10;
    const n: usize = 1024;

    // decode encoded signature
    let r = dec_sig(logn, &signature);

    let salt = r.0;
    let s1 = r.1;

    // decode public key
    let r = dec_pub(logn, &pub_key);
    let q00 = r.0;
    let q01 = r.1;

    // check if public key decoding failed
    if q00[0] == 0 && q00.len() == 1 && q01[0] == 0 && q01.len() == 1 {
        println!("failed on decoding public key");
        return false;
    }

    // compute hash of public key
    let mut hpubshaker = Shake256::default();
    hpubshaker.update(&pub_key[..]);
    let mut hpub: [u8; LENHPUB] = [0; LENHPUB];
    hpubshaker.finalize_xof_reset_into(&mut hpub);

    // compute hash M
    let mut shaker = Shake256::default();
    shaker.update(msg);
    shaker.update(&hpub);
    let mut m: [u8; 64] = [0; 64];
    shaker.finalize_xof_reset_into(&mut m);

    shaker.update(&m);
    shaker.update(&salt);

    let mut h: [u8; n / 4] = [0; n / 4];
    shaker.finalize_xof_reset_into(&mut h);

    // convert h to two polynomials
    let (h0, h1) = (
        &bytes_to_poly(&h[0..n / 8], n),
        &bytes_to_poly(&h[(n / 8)..n / 4], n),
    );

    let w1 = poly_sub(&h1, &poly_times_const(&i16vec_to_i32vec(&s1), 2));

    if !symbreak(&w1) {
        println!("Symbreak failed");
        return false;
    }

    let q00_i32 = i16vec_to_i32vec(&q00);
    let q01_i32 = i16vec_to_i32vec(&q01);
    let q00_i64 = i32vec_to_i64vec(&q00_i32);
    let q01_i64 = i32vec_to_i64vec(&q01_i32);
    let w1_i32 = i64vec_to_i32vec(&w1);
    let h0_i32 = i64vec_to_i32vec(&h0);

    let w0 = rebuildw0(logn, &q00_i32, &q01_i32, &w1_i32, &h0_i32);

    if w0[0] == 0 && w0.len() == 1 {
        println!("failed on rebuild");
        return false;
    }

    let (p1, p2): (i64, i64) = (2147473409, 2147389441);

    let r1 = poly_qnorm(&q00_i64, &q01_i64, &i32vec_to_i64vec(&w0), &w1, p1);
    let r2 = poly_qnorm(&q00_i64, &q01_i64, &i32vec_to_i64vec(&w0), &w1, p2);

    if r1 != r2 || modulo(r1, n as i64) != 0 {
        println!("failed here");
        return false;
    }

    let r1 = r1 / (n as i64);

    let sigmaverify: f64 = 1.042;

    if (r1 as f64) > (8 * n) as f64 * sigmaverify.powi(2) {
        return false;
    }

    return true;
}

fn i64vec_to_i32vec(f: &Vec<i64>) -> Vec<i32> {
    let res: Vec<i32> = f.iter().map(|&x| x as i32).collect();

    return res;
}

fn i16vec_to_i32vec(f: &Vec<i16>) -> Vec<i32> {
    let res: Vec<i32> = f.iter().map(|&x| x as i32).collect();
    return res;
}

fn i32vec_to_i64vec(f: &Vec<i32>) -> Vec<i64> {
    let res: Vec<i64> = f.iter().map(|&x| x as i64).collect();
    return res;
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
