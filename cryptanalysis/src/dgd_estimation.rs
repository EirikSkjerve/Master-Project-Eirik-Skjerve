// this file contains experiments for estimating the std. dev. and variance
// for the practical implementation of the "Discrete Gaussian Distribution" in Hawk.

use hawklib::hawkkeygen::hawkkeygen;
use hawklib::hawksign::hawksign_total;

use rand::Rng;
use nalgebra::*;

fn get_random_bytes(num_bytes: usize) -> Vec<u8> {
    // return num_bytes random bytes in a Vec<u8>
    //
    // example: get_random_bytes(4) -> vec![100,200,33,99]

    let mut res = Vec::with_capacity(num_bytes);
    let mut rng = rand::thread_rng();

    for _ in 0..num_bytes {
        res.push(rng.gen_range(0..255));
    }

    res
}

pub fn estimate_sigma(t: usize, n: usize) {
    // create t x-vectors with hawk degree n
    
    assert!(n==256 || n==512 || n==1024);

    // generate a keypair
    let (privkey, pubkey) = hawkkeygen(n);

    // create t messages
    let mut messages: Vec<Vec<u8>> = Vec::with_capacity(t);
    println!("Generating {t} messages...");
    for _ in 0..t {
        messages.push(get_random_bytes(100));
    }

    // create collection of t samples corresponding to the above messages
    println!("Generating {t} samples...");
    let mut xsamples: Vec<Vec<i64>> = Vec::with_capacity(t);
    for i in 0..t {
        // hawksign_total returns the raw x-vector as the second element in tuple
        xsamples.push(hawksign_total(&privkey, &messages[i], n).1);
    }

    println!("{}, {}", xsamples.len(), xsamples[0].len());

    let xsamples_flat: Vec<i64> = xsamples.into_iter().flatten().collect();
    println!("{}", xsamples_flat.len());
}
