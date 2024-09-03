use sha3::{
    digest::{ExtendableOutput, ExtendableOutputReset, Update},
    Shake256,
};

use crate::fft::{div_fft_f64, mul_fft_f64};
use crate::ntt::*;
use crate::sign::symbreak;
use crate::utils::{adjoint, bytes_to_poly, mod_pow, modulo, poly_add, poly_mult_ntt, poly_sub};
use crate::verifyutils::*;

pub fn verify(
    logn: usize,
    msg: usize,
    q00: Vec<i64>,
    q01: Vec<i64>,
    signature: (Vec<u8>, Vec<i64>)
) -> bool {
    let n = 1 << logn;

    let salt = signature.0;
    let s1 = signature.1;

    // compute hash M
    let mut shaker = Shake256::default();
    shaker.update(&msg.to_ne_bytes());
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

    println!("h0: {:?} \nh1: {:?}", h0, h1);

    let w1 = poly_sub(&h1, &poly_times_const(&s1, 2));

    if !symbreak(&w1) {
        println!("Symbreak failed");
        return false;
    }


    let q00_i32 = i64vec_to_i32vec(&q00);
    let q01_i32 = i64vec_to_i32vec(&q01);
    let w1_i32 = i64vec_to_i32vec(&w1);
    let h0_i32 = i64vec_to_i32vec(&h0);
    let w0 = rebuildw0(logn, &q00_i32, &q01_i32, &w1_i32, &h0_i32);

    let (p1, p2): (i64, i64) = (2147473409, 2147389441);
    // println!("q00 = {:?} \nq01 = {:?} \nw0 = {:?} \nw1 = {:?}", q00, q01, w0, w1);
    let r1 = polyQnorm(&q00, &q01, &i32vec_to_i64vec(&w0), &w1, p1);
    let r2 = polyQnorm(&q00, &q01, &i32vec_to_i64vec(&w0), &w1, p2);

    println!("r1: {} \nr2: {}", r1, r2);
    if r1 != r2 || modulo(r1, n as i64) != 0 {
        // println!("failed here");
        return false;
    }

    let r1 = r1 / n as i64;

    let sigmaverify: f64 = 1.042;

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

fn i32vec_to_i64vec(f: &Vec<i32>) -> Vec<i64> {
    let res: Vec<i64> = f.iter().map(|&x| x as i64).collect();
    return res;
}

pub fn vec_to_slice(vec: &Vec<u8>) -> &[u8] {
    vec
}
//
// pub fn rebuilds0(
//     q00: &Vec<i64>,
//     q01: &Vec<i64>,
//     w1: &Vec<i64>,
//     h0: &Vec<i64>,
//     h1: &Vec<i64>,
// ) -> Vec<i64> {
//     let a = poly_div_const_f64(h0, 2);
//     let b = div_fft_f64(q01, q00);
//     let c = poly_div_const_f64(h1, 2);
//     let w1_f64: Vec<f64> = w1.iter().map(|&x| x as f64).collect();
//     let d = poly_sub_f64(&c, &w1_f64);
//
//     let s0_f64 = poly_add_f64(&a, &mul_fft_f64(&b, &d));
//
//     // round each coefficient
//     let s0 = s0_f64.iter().map(|&x| x.round() as i64).collect();
//
//     // compute w0 = h0 - 2*s0
//     let w0 = poly_sub(&h0, &poly_times_const(&s0, 2));
//     return w0;
// }
//

