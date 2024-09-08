use sha3::{
    digest::{ExtendableOutput, ExtendableOutputReset, Update},
    Shake256,
};

use crate::codec::{dec_pub, dec_sig};
use crate::sign::symbreak;
use crate::utils::{bytes_to_poly, modulo, poly_sub};
use crate::verifyutils::*;

pub fn verify(logn: usize, msg: &str, pub_key: &Vec<u8>, signature: &Vec<u8>) -> bool {
    let n = 1 << logn;

    // decode encoded signature
    let r = dec_sig(logn, &signature);

    let salt = r.0;
    let s1 = r.1;

    // println!("sig from verify: {:?}", s1);

    // println!("salt from verify: {:?}", salt);
    // println!("s1 from verify: {:?}", s1);
    // TODO check signature

    // decode public key
    let r = dec_pub(logn, &pub_key);
    let q00 = r.0;
    let q01 = r.1;
    // println!("q00 = {:?}", q00);
    // println!("q01 = {:?}", q01);

    // check if public key decoding failed
    if q00[0] == 0 && q00.len() == 1 && q01[0] == 0 && q01.len() == 1 {
        println!("failed on decoding public key");
        return false;
    }

    // compute hash of public key
    let mut hpubshaker = Shake256::default();
    hpubshaker.update(&pub_key[..]);
    let mut hpub: [u8; 16] = [0; 16];
    hpubshaker.finalize_xof_reset_into(&mut hpub);

    // compute hash M
    let mut shaker = Shake256::default();
    shaker.update(&msg.as_bytes());
    shaker.update(&hpub);
    let mut m: [u8; 64] = [0; 64];
    shaker.finalize_xof_reset_into(&mut m);

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


    let w1 = poly_sub(&h1, &poly_times_const(&i16vec_to_i32vec(&s1), 2));
    // println!("w1 <- h1 - 2*s1: \n {:?} <- \n{:?} -")

    println!("w1 from verify: {:?}", w1);

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
    // println!("q00 = {:?}", q00_i32);
    // println!("q01 = {:?}", q01_i32);
    // println!("w1 = {:?}", w1_i32);
    // println!("h0 = {:?}", h0_i32);
    //
    // println!("\n w0 = {:?}", w0);
    if w0[0] == 0 && w0.len() == 1 {
        println!("failed on rebuild");
        return false;
    }
    // println!("rebuilt w0: {:?}", w0);

    let (p1, p2): (i64, i64) = (2147473409, 2147389441);
    // println!("q00 = {:?} \nq01 = {:?} \nw0 = {:?} \nw1 = {:?}", q00, q01, w0, w1);
    // println!("h0 = {:?}", h0);

    let r1 = polyQnorm(&q00_i64, &q01_i64, &i32vec_to_i64vec(&w0), &w1, p1);
    let r2 = polyQnorm(&q00_i64, &q01_i64, &i32vec_to_i64vec(&w0), &w1, p2);

    // println!("q00: {:?} \nq01: {:?}", q00_i32, q01_i32);
    // println!("r1: {} \nr2: {}", r1, r2);

    if r1 != r2 || modulo(r1, n as i64) != 0 {
        println!("failed here");
        return false;
    }

    let r1 = r1 / (n as i64);

    let sigmaverify: f64 = 1.042;

    println!(
        "r1: {} \n limit: {}",
        r1,
        (8 * n) as f64 * sigmaverify.powi(2)
    );
    if (r1 as f64) > (8 * n) as f64 * sigmaverify.powi(2) {
        println!("Too big");
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

pub fn vec_to_slice(vec: &Vec<u8>) -> &[u8] {
    vec
}
