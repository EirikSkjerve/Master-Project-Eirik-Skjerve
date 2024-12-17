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

pub fn estimate_sigma_time(t: usize, n: usize) {
    // create t x-vectors with hawk degree n
    // this function will keep the t samples in memory and reuse them for the estimation of mu and
    // sigma
    // Drawback is that there will be an upper limit to how many samples can be kept in memory at
    // given time.

    assert!(n == 256 || n == 512 || n == 1024);

    // generate a keypair
    let (privkey, _) = hawkkeygen(n);

    println!("Estimating sigma sampling {t} vectors of length 2*{n}");

    let start = Instant::now();

    let mut samples: Vec<Vec<i64>> = Vec::new();
    println!("Generating samples...");
    // estimating mu
    for _ in 0..t {
        // sample x-vector of length 2n given a random vector
        samples.push(hawksign_x_only(&privkey, &get_random_bytes(100), n));
    }

    let mut end = start.elapsed();
    println!("\nTime used: {:?}\n", end);

    let start = Instant::now();
    println!("Estimating mu...");
    // compute an estimate mu
    let mut mu: f64 = samples.iter().flatten().map(|&x| x as f64).sum();
    mu /= (t * 2 * n) as f64;

    end = start.elapsed();
    println!("Time used: {:?}\n", end);

    let start = Instant::now();
    println!("Estimating var...");
    // compute an estimate of sigma^2 using estimated mu
    let mut var: f64 = samples
        .iter()
        .flatten()
        .map(|&x| (x as f64 - mu).powi(2))
        .sum();

    var /= (t * 2 * n - 1) as f64;

    end = start.elapsed();
    println!("Time used: {:?}\n", end);

    println!("\nResults of estimating:");
    println!("Exp[X]: {mu}");
    println!("Var[X]: {var}");
    println!("Sigma: {}", var.sqrt());

}


pub fn estimate_sigma_mem(t: usize, n: usize) {
    // create t x-vectors with hawk degree n
    // this function will not store any intermediate samples and therefore will not need a lot of
    // memory. Drawback is that t x-vectors have to be generated twice; once for mu and once for
    // var
    //
    // TODO need to check that rng seeds are properly set up with respect to number of samples we
    // want to collect

    assert!(n == 256 || n == 512 || n == 1024);

    // generate a keypair
    let (privkey, _) = hawkkeygen(n);

    println!("Estimating sigma sampling {t} vectors of length 2*{n}");

    let start = Instant::now();

    println!("Generating samples for mu...");
    // estimating mu
    let mut mu: f64 = 0.0;
    for _ in 0..t {
        // sample x-vector of length 2n given a random vector
        let temp: i64 = hawksign_x_only(&privkey, &get_random_bytes(100), n)
            .iter()
            .sum();
        mu += temp as f64 /(t*2*n) as f64;
    }

    // mu /= (t * 2 * n) as f64;
    println!("Time used estimating mu: {:?}\n", start.elapsed());

    let start = Instant::now();
    println!("Generating samples for var...");
    // estimating var
    let mut var: f64 = 0.0;
    for _ in 0..t {
        // sample x-vector of length 2n given a random vector
        let temp: f64 = hawksign_x_only(&privkey, &get_random_bytes(100), n)
            .iter()
            .map(|&x| (x as f64 - mu).powi(2))
            .sum();
        var += temp / (t*2*n) as f64;
    }

    // var /= (t * 2 * n) as f64;
    println!("Time used estimating var: {:?}\n", start.elapsed());


    println!("\nResults of estimating:");
    println!("Exp[X]:  {mu}");
    println!("Var[X]:  {var}");
    println!("Sigma:   {}", var.sqrt());
    println!("Sigma/2: {}", var.sqrt()/2.0);

}
