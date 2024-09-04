use crate::compress::compressgr;
use crate::decompress::decompressgr;
use crate::grutils::*;
use crate::params::{params_f, params_i};
use crate::utils::bytes_to_poly;

pub fn enc_pub(logn: usize, q00: &Vec<i64>, q01: &Vec<i64>) -> Vec<u8> {
    /*
     * encode public key
     * failure to encode yields [0] as a result
     */

    let n = 1 << logn;
    if q00[0] < -(1 << 15) || q00[0] >= (1 << 15) {
        println!("failure from encpub 1");
        // default failure value
        return vec![0];
    }

    let v: usize = 16 - params_i(logn, "high00") as usize;

    let mut q00_c = q00.clone();
    q00_c[0] = (q00[0]) >> v;

    let mut y00 = compressgr(
        &q00_c[0..(n / 2)].to_vec(),
        params_i(logn, "low00") as usize,
        params_i(logn, "high00") as usize,
    );

    if y00[0] == 0 && y00.len() == 1 {
        println!("failure from encpub 2");
        // this is the failure return value
        return vec![0];
    }

    let bin_temp = bin(modulo(q00[0], 1 << v) as i32, v);
    for b_t in bin_temp {
        y00.push(b_t);
    }
    while modulo(y00.len(), 8) != 0 {
        y00.push(0);
    }

    let y01 = compressgr(
        q01,
        params_i(logn, "low01") as usize,
        params_i(logn, "high01") as usize,
    );

    if y01[0] == 0 && y01.len() == 1 {
        // failure return value
        println!("failure from encpub 3");
        return vec![0];
    }

    let mut y: Vec<u8> = Vec::with_capacity(y00.len() + y01.len());
    for i in 0..y00.len() {
        y.push(y00[i]);
    }
    for i in 0..y01.len() {
        y.push(y01[i]);
    }

    if y.len() > (params_i(logn, "lenpub") * 8) as usize {
        // failure return value
        println!("failure from encpub 4");
        return vec![0];
    }
    println!("y len: {}", y.len());

    // padding of the y vector
    while y.len() < (params_i(logn, "lenpub") * 8) as usize {
        y.push(0);
    }
    println!("y len: {}", y.len());
    let packed = packbits(&y);

    return packed;
}

pub fn dec_pub(logn: usize, pub_enc: &Vec<u8>) -> (Vec<i16>, Vec<i16>) {
    /*
     * Decodes an encoded public key
     * Serves as the inverse function of enc_pub()
     * failure returns ([0], [0])
     */

    let n = 1 << logn;

    if pub_enc.len() != params_i(logn, "lenpub") as usize {
        // failure return value
        println!("failure from decpub 1");
        return (vec![0], vec![0]);
    }

    let v = 16 - params_i(logn, "high00") as usize;
    let y = unpackbits(&pub_enc);

    let r00 = decompressgr(
        &y,
        n / 2,
        params_i(logn, "low00") as usize,
        params_i(logn, "high00") as usize,
    );
    if r00.0[0] == 0 && (r00.0).len() == 1 {
        println!("failure from decpub 2");
        return (vec![0], vec![0]);
    }

    let (r00, j) = (r00.0, r00.1);

    let mut q00: Vec<i16> = vec![0; n];

    for i in 0..(n / 2) {
        q00[i] = r00[i] as i16;
    }

    if y.len() * 8 < j + v {
        println!("failure from decpub 3");
        return (vec![0], vec![0]);
    }

    q00[0] *= 1 << v;
    q00[0] += int(&y[j..(j + v)].to_vec()) as i16;

    let mut j = j + v;

    while modulo(j, 8) != 0 {
        if j >= y.len() || y[j] != 0 {
            println!("failure from decpub 4");
            return (vec![0], vec![0]);
        }
        j += 1;
    }

    q00[n / 2] = 0;
    for i in ((n / 2) + 1)..n {
        q00[i] = -q00[n - i];
    }

    let r01 = decompressgr(
        &y[j..y.len()].to_vec(),
        n,
        params_i(logn, "low01") as usize,
        params_i(logn, "high01") as usize,
    );

    if r01.0[0] == 0 && r01.0.len() == 1 {
        println!("failure from decpub 5");
        return (vec![0], vec![0]);
    }

    let (r01, jp) = (r01.0, r01.1);

    j += jp;
    let mut q01 = Vec::with_capacity(r01.len());
    for i in 0..r01.len() {
        q01.push(r01[i] as i16);
    }

    while j < y.len() {
        if y[j] != 0 {
            println!("failure from decpub 6");
            return (vec![0], vec![0]);
        }
        j += 1;
    }

    return (q00, q01);
}

pub fn enc_priv() {}

pub fn dec_priv() {}

pub fn enc_sig() {}

pub fn dec_sig() {}
