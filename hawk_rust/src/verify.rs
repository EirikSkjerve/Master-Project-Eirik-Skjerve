use sha3::{
    digest::{ExtendableOutput, ExtendableOutputReset, Update},
    Shake256,
};

use crate::sign::symbreak;
use crate::utils::{bytes_to_poly, modulo, poly_add, poly_mult_ntt, poly_sub};
use crate::fft::{div_fft_f64, mul_fft_f64};
use crate::ntt::*;

pub fn verify(msg: usize, q00: Vec<i64>, q01: Vec<i64>, signature: (Vec<u8>, Vec<i64>), logn: usize) {

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
    let mut h: [u8; 256/4] = [0; 256/4];
    shaker.finalize_xof_reset_into(&mut h);

    // convert h to two polynomials
    let (h0, h1) = (&bytes_to_poly(&h[0..256/8], n), &bytes_to_poly(&h[(256/8)..256/4], n));

    let w1 = poly_sub(&h1, &poly_times_const(&s1, 2));

    if !symbreak(&w1) {
        return;
    }

    let w0 = rebuilds0(&q00, &q01, &w1, &h0, &h1);

    let (p1, p2): (i64, i64) = (2147473409, 2147389441);

 
}

pub fn vec_to_slice(vec: &Vec<u8>) -> &[u8] {
    vec
}
pub fn poly_div_const_i64(f: &Vec<i64>, c: i64) -> Vec<i64> {
    return f.iter().map(|&x| x/c).collect();
}

fn rebuilds0(q00: &Vec<i64>, q01: &Vec<i64>, w1: &Vec<i64>, h0: &Vec<i64>, h1: &Vec<i64>) -> Vec<i64> {


    let a = poly_div_const_f64(h0, 2);
    let b = div_fft_f64(q01, q00);
    let c = poly_div_const_f64(h1, 2);
    let w1_f64: Vec<f64> = w1.iter().map(|&x| x as f64).collect();
    let d = poly_sub_f64(&c, &w1_f64);

    let s0_f64 = poly_add_f64(&a, &mul_fft_f64(&b, &d));

    // round each coefficient
    let s0 = s0_f64.iter().map(|&x| x.round() as i64).collect();

    return s0;
}

fn polyQnorm(q00: &Vec<i64>, q01: &Vec<i64>, w0: &Vec<i64>, w1: &Vec<i64>, p: i64) {

}

pub fn poly_times_const(f: &Vec<i64>, c: i64) -> Vec<i64> {
    return f.iter().map(|&x| x*c).collect();
}

pub fn poly_div_const_f64(f: &Vec<i64>, c: i64) -> Vec<f64> {
    return f.iter().map(|&x| x as f64 / c as f64).collect();
}

fn poly_sub_f64(f: &Vec<f64>, g: &Vec<f64>) -> Vec<f64> {
    let mut h: Vec<f64> = vec![0.0; f.len()];
    for i in 0..h.len(){
        h[i] = f[i] - g[i];
    }
    return h;
}

fn poly_add_f64(f: &Vec<f64>, g: &Vec<f64>) -> Vec<f64> {
    let mut h: Vec<f64> = vec![0.0; f.len()];
    for i in 0..h.len(){
        h[i] = f[i] + g[i];
    }
    return h;
}
