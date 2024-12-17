// this file contains experiments for estimating the std. dev. and variance
// for the practical implementation of the "Discrete Gaussian Distribution" in Hawk.

use hawklib::hawkkeygen::hawkkeygen;
use hawklib::hawksign::{hawksign_total, hawksign_x_only};

use rand::Rng;
use std::time::Instant;

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

    assert!(n == 256 || n == 512 || n == 1024);

    // generate a keypair
    let (privkey, _) = hawkkeygen(n);

    println!("Estimating sigma sampling {t} vectors of length 2*{n}");

    let start = Instant::now();

    let mut samples: Vec<Vec<i64>> = Vec::new();

    // estimating mu
    for _ in 0..t {
        // sample x-vector of length 2n given a random vector
        samples.push(hawksign_x_only(&privkey, &get_random_bytes(100), n));
    }

    // compute an estimate mu
    let mut mu: f64 = samples.iter().flatten().map(|&x| x as f64).sum();

    mu /= (t * 2 * n) as f64;

    // compute an estimate of sigma^2 using estimated mu
    let mut var: f64 = samples
        .iter()
        .flatten()
        .map(|&x| (x as f64 - mu).powi(2))
        .sum();

    var /= (t * 2 * n - 1) as f64;

    let end = start.elapsed();

    println!("\nResults of estimating:");
    println!("Exp[X]: {mu}");
    println!("Var[X]: {var}");
    println!("Sigma: {}", var.sqrt());

    println!("\nTime used: {:?}", end);
}
