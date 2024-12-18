// this file contains experiments for estimating the std. dev. and variance
// for the practical implementation of the "Discrete Gaussian Distribution" in Hawk.

use hawklib::hawkkeygen::hawkkeygen;
use hawklib::hawksign::{hawksign_total, hawksign_x_only};
use hawklib::parameters::{hawk1024_params, hawk256_params, hawk512_params};

use rand::Rng;
use std::time::{Instant, Duration};
use std::fs::File;

use prettytable::{Cell, Row, Table, Attr, color};

pub fn estimate_sigma_mem_all(t: usize, store_file: bool){
    let ns = vec![256, 512, 1024];

    let precision = 5;

    // define table of estimations
    let mut table = Table::new();

    table.add_row(row![i->"deg", i->"Mu", i->"Sigma^2", i->"Sigma", i->"Th. Sigma", i->"Time (s)"]);

    for n in ns {

        // get the right parameters. These are just for comparison
        let sigmasign = match n {
            256 => hawk256_params::SIGMASIGN,
            512 => hawk512_params::SIGMASIGN,
            _ => hawk1024_params::SIGMASIGN,
        };

        // run the estimation
        let (mu, sig, time) = estimate_sigma_mem(t, n);

        // write results to table
        table.add_row(
            row![
            FG->n.to_string(),
            FM->format!("{:.1$}", mu, precision), 
            Fr->format!("{:.1$}", sig, precision),
            Fw->format!("{:.1$}", sig.sqrt(), precision),
            Fc->format!("{}", 2.0*sigmasign),
            Fy->format!("{:.1$}", time.as_secs_f64().to_string(), precision)
            ]);
    }
    table.printstd();

    if store_file {
        let out = File::create("table_csv.csv").unwrap();
        table.to_csv(out).unwrap();
        println!("Created file at table_csv.csv");
    }

}

pub fn estimate_sigma_mem(t: usize, n: usize) -> (f64, f64, Duration){
    // create t x-vectors with hawk degree n
    // this function will not store any intermediate samples and therefore will not need a lot of
    // memory. Drawback is that t x-vectors have to be generated twice; once for mu and once for
    // var
    //
    // returns (Exp[x], Var[X], time used)
    //


    assert!(n == 256 || n == 512 || n == 1024);


    let no_retry = false;

    // generate a keypair
    let (privkey, _) = hawkkeygen(n);

    println!("Estimating sigma sampling {t} vectors of length 2*{n}");

    // measure time for mu
    let start = Instant::now();

    // estimating mu
    let mut mu: f64 = 0.0;
    for _ in 0..t {
        // sample x-vector of length 2n given a random vector
        let temp: i64 = hawksign_x_only(&privkey, &get_random_bytes(100), n, no_retry)
            .iter()
            .sum();
        mu += temp as f64 /(t*2*n) as f64;
    }

    // estimating var
    let mut var: f64 = 0.0;
    for _ in 0..t {
        // sample x-vector of length 2n given a random vector
        let temp: f64 = hawksign_x_only(&privkey, &get_random_bytes(100), n, no_retry)
            .iter()
            .map(|&x| (x as f64 - mu).powi(2))
            .sum();
        var += temp / (t*2*n) as f64;
    }

    let end = start.elapsed();

    (mu, var, end)

    // return data on the following form:

}

pub fn estimate_sigma_time(t: usize, n: usize) {
    // create t x-vectors with hawk degree n
    // this function will keep the t samples in memory and reuse them for the estimation of mu and
    // sigma
    // Drawback is that there will be an upper limit to how many samples can be kept in memory at
    // given time.

    let no_retry = false;

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
        samples.push(hawksign_x_only(&privkey, &get_random_bytes(100), n, no_retry));
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
