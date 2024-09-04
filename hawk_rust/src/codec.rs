use crate::utils::bytes_to_poly;
use crate::params::{params_f, params_i};
use crate::compress::compressgr;
use crate::decompress::decompressgr;
use crate::grutils::*;

pub fn enc_pub(logn: usize, q00: &Vec<i64>, q01: &Vec<i64>) -> Vec<u8> {

    let n = 1 << logn;
    if q00[0] < -(1<<15) || q00[0] >= (1<<15){
        // default failure value
        return vec![0]
    }

    let v: usize = 16 - params_i(logn, "high00") as usize;
    let mut q00_c = q00.clone();
    q00_c[0] = (q00[0])>>v;
    
    let mut y00 = compressgr(&q00_c[0..(n/2)].to_vec(), params_i(logn, "low00") as usize, params_i(logn, "high00") as usize);

    if y00[0] == 0 && y00.len() == 1 {
        // this is the failure return value
        return vec![0];
    }

    let bin_temp = bin(modulo(q00[0], 1<<v) as i32, v);
    for b_t in bin_temp{
        y00.push(b_t);
    }

    let y01 = compressgr(q01, params_i(logn, "low01") as usize, params_i(logn, "high01") as usize);

    if y01[0] == 0 && y00.len() == 1 {
        // failure return value
        return vec![0];
    }
    
    let mut y: Vec<u8> = Vec::with_capacity(y00.len() + y01.len());
    for i in 0..y00.len(){
        y.push(y00[i]);
    }

    for i in 0..y01.len(){
        y.push(y01[i]);
    }

    if y.len() > (params_i(logn, "lenpub")*8) as usize {
        // failure return value
        return vec![0];
    }

    // padding of the y vector
    while y.len() < (params_i(logn, "lenpub")*8) as usize {
        y.push(0);
    }

    let packed = packbits(&y);

    return packed;
}

pub fn dec_pub() {

}

pub fn enc_priv() {

}

pub fn dec_priv() {

}

pub fn enc_sig() {

}

pub fn dec_sig() {

}
