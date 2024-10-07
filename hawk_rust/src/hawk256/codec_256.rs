// CODEC for hawk 256

use crate::compression::compress::compressgr;
use crate::compression::decompress::decompressgr;
use crate::grutils::*;
use crate::utils::bytes_to_poly;

use crate::parameters::hawk256_params::*;

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

    let v: usize = 16 - HIGH00 as usize;

    let mut q00_c = q00.clone();
    q00_c[0] = (q00[0]) >> v;

    let mut y00 = compressgr(
        &q00_c[0..(n / 2)].to_vec(),
        LOW00,
        HIGH00,
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
        LOW01,
        HIGH01,
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

    if y.len() > (LENPUB * 8) as usize {
        // failure return value
        println!("failure from encpub 4");
        return vec![0];
    }

    // padding of the y vector
    while y.len() < (LENPUB * 8) as usize {
        y.push(0);
    }
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

    if pub_enc.len() != LENPUB {
        // failure return value
        println!("failure from decpub 1");
        return (vec![0], vec![0]);
    }

    let v = 16 - HIGH00;
    let y = unpackbits(&pub_enc);

    let r00 = decompressgr(
        &y,
        n / 2,
        LOW00,
        HIGH00,
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
        LOW01,
        HIGH01,
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

pub fn enc_priv(kgseed: &[u8], bigf: &Vec<i64>, bigg: &Vec<i64>, hpub: &[u8]) -> Vec<u8> {
    let mut bigfbiggmod2 = Vec::with_capacity(bigf.len() + bigg.len());

    // compute bigf mod 2 and bigg mod 2 and append them into a single vector
    for i in 0..bigf.len() {
        bigfbiggmod2.push(modulo(bigf[i], 2) as u8);
    }
    for i in 0..bigg.len() {
        bigfbiggmod2.push(modulo(bigg[i], 2) as u8);
    }

    // initialize result vector
    let mut res: Vec<u8> = Vec::new();

    // append bytes of kgseed to result vector
    for kb in kgseed.iter() {
        res.push(*kb);
    }

    // pack the bits in the joined bigfmod2||biggmod2
    let packedfg = packbits(&bigfbiggmod2);

    // add the packed bits into result vector
    for pfg in packedfg.iter() {
        res.push(*pfg);
    }

    // add the bytes from hpub
    for hp in hpub.iter() {
        res.push(*hp);
    }

    // return the result
    return res;
}

pub fn dec_priv(logn: usize, priv_enc: &Vec<u8>) -> (Vec<u8>, Vec<i64>, Vec<i64>, Vec<u8>) {
    let n = 1 << logn;
    // make input encoded private key to a vector
    let priv_vec = priv_enc.to_vec();
    // get the length of kgseed
    let lenkgseed = LENKGSEED;
    let kgseed_arr = &priv_vec[0..lenkgseed];

    // vectors for the polynomials bigf mod 2 and bigg mod 2
    let bigfmod2 = bytes_to_poly(&priv_vec[lenkgseed..lenkgseed + (n / 8)], 1 << logn);
    let biggmod2 = bytes_to_poly(
        &priv_vec[lenkgseed + (n / 8)..lenkgseed + (n / 4)],
        1 << logn,
    );

    let lenhpub = LENHPUB;
    let hpub = &priv_vec[(priv_vec.len() - lenhpub)..priv_vec.len()].to_vec();

    return (
        kgseed_arr.to_vec(),
        bigfmod2.to_vec(),
        biggmod2.to_vec(),
        hpub.clone(),
    );
}

pub fn enc_sig(logn: usize, salt: &Vec<u8>, s1: &Vec<i64>) -> Vec<u8> {
    // compress s1
    let mut y = compressgr(
        &s1,
        LOWS1,
        HIGHS1,
    );

    // check if compression has failed
    if y[0] == 0 && y.len() == 1 {
        println!("failure from enc sig 1");
        return vec![0];
    }

    // get the prescribed length of y
    let leny = ((LENSIG - LENSALT) * 8) as usize;

    // return failure if y is longer than prescribed
    if y.len() > leny {
        println!("failure from enc sig 2");
        return vec![0];
    }

    // pad y with 0's
    while y.len() < leny {
        y.push(0);
    }

    // append salt || packedbits(y)
    let mut res = salt.clone();
    let y_packed = packbits(&y);
    for yp in y_packed.iter() {
        res.push(*yp);
    }

    return res;
}

pub fn dec_sig(logn: usize, sig_enc: &Vec<u8>) -> (Vec<u8>, Vec<i16>) {
    let n = 1 << logn;

    // check length of encoded signature
    if sig_enc.len() != LENSIG {
        return (vec![0], vec![0]);
    }

    let salt = &sig_enc[0..LENSALT];

    // unpack encoded sig to bits
    let y = unpackbits(&sig_enc);

    let s1 = decompressgr(
        &y[((LENSALT * 8) as usize)..y.len()].to_vec(),
        n,
        LOWS1,
        HIGHS1,
    );

    // unpack the values in the tuple
    let mut j = s1.1;
    let s1 = s1.0;

    // return failure value if decompress fails
    if s1[0] == 0 && s1.len() == 1 {
        return (vec![0], vec![0]);
    }

    let s1: Vec<i16> = s1.iter().map(|&x| x as i16).collect();

    j += (LENSALT * 8) as usize;

    while j < y.len() {
        if y[j] != 0 {
            return (vec![0], vec![0]);
        }
        j += 1;
    }

    return (salt.to_vec(), s1);
}
