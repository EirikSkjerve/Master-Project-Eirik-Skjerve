// this file contains experiments for estimating the std. dev. and variance
// for the practical implementation of the "Discrete Gaussian Distribution" in Hawk.

use hawklib::hawkkeygen::hawkkeygen;
use hawklib::hawksign::{hawksign_total, hawksign_x_only};
use hawklib::parameters::{hawk1024_params, hawk256_params, hawk512_params};

use rand::Rng;
use std::fs::File;
use std::path::Path;
use std::time::{Duration, Instant};
use std::io::{stdout, Write};

use prettytable::{color, Attr, Cell, Row, Table};

pub fn estimate_sigma_mem_all(t: usize, store_file: bool) {
    let ns = vec![256, 512, 1024];
    // let ns = vec![256];

    let precision = 8;

    // define table of estimations
    let mut table = Table::new();

    table.add_row(row![
                  i->"Deg",
                  i-> "Num\nVectors",
                  i->"Mu",
                  i->"Var",
                  i->"Sigma",
                  i->"Th. \nSigma",
                  i->"Diff.",
                  i->"Kur",
                  i->"3 Var ^2",
                  i->"Diff.",
                  i->"Time (s)",
    ]);

    for n in ns {
        // get the right parameters. These are just for comparison
        let sigmasign = match n {
            256 => hawk256_params::SIGMASIGN,
            512 => hawk512_params::SIGMASIGN,
            _ => hawk1024_params::SIGMASIGN,
        };

        // run the estimation
        let (mu, var, kur, time) = estimate_sigma_mem(t, n);

        // write results to table
        table.add_row(row![
        FG->n.to_string(),
        FW->(t).to_string(),
        FM->format!("{:.1$}", mu, precision),
        FR->format!("{:.1$}", var, precision),
        Fw->format!("{:.1$}", var.sqrt(), precision),
        Fc->format!("{}", 2.0*sigmasign),
        Fy->format!("{:.1$}", (var.sqrt() - 2.0*sigmasign).abs(), precision),
        FR->format!("{:.1$}", kur, precision),
        Fc->format!("{:.1$}", 3.0*var.powi(2), precision),
        Fy->format!("{:.1$}", (kur- 3.0*var.powi(2)).abs(), precision),
        Fw->format!("{:.1$}", time.as_secs_f64().to_string(), precision)
        ]);
    }
    table.printstd();

    // store file
    if store_file {
        // making sure a file is never written over
        let pathname_base = String::from("sigma_table");
        let mut pathname = String::from("sigma_table_0");
        let mut ctr = 1;

        loop {
            if Path::new(&format!("{}.csv", pathname)).is_file() {
                // println!("{}.csv already exists!", pathname);
                pathname = format!("{}_{}", pathname_base, ctr);
                ctr += 1;
                continue;
            }
            break;
        }

        let out = File::create(&format!("{}.csv", pathname)).unwrap();
        table.to_csv(out).unwrap();
        println!("Created file at {}.csv", pathname);
    }
}

pub fn estimate_sigma_mem(t: usize, n: usize) -> (f64, f64, f64, Duration) {
    // create t x-vectors with hawk degree n
    // this function will not store any intermediate samples and therefore will not need a lot of
    // memory. Drawback is that t x-vectors have to be generated twice; once for mu and once for
    // var
    //
    // returns (Exp[x], Var[X], time used)
    //

    // generate a keypair
    let (privkey, _) = hawkkeygen(n);

    assert!(n == 256 || n == 512 || n == 1024);

    // for nice printouts
    let mut stdout = stdout();

    // bool determining if sampling of x should retry if a sampled x is not "valid"
    // this has an effect on the measurement of the practical distribution
    let no_retry = false;

    println!(
        "\nEstimating sigma by sampling {} vectors of length 2*{n}",
        t
    );

    // measure time for mu
    let start = Instant::now();

    // estimating mu
    let mut mu: f64 = 0.0;
    for i in 0..t {

        // stdout.flush().unwrap();
        // print!("\r{} %", (((i+1) as f64 / t as f64) * 100.0).abs());

        // sample x-vector of length 2n given a random vector
        let temp: i64 = hawksign_x_only(&privkey, &get_random_bytes(100), n, no_retry)
            .iter()
            .sum();
        mu += temp as f64 / (t * 2 * n) as f64;
    }

    // estimating var
    let mut var: f64 = 0.0;
    let mut kur: f64 = 0.0;

    for i in 0..t {

        // stdout.flush().unwrap();
        // print!();

        // sample x-vector of length 2n given a random vector
        let temp: Vec<i64> = hawksign_x_only(&privkey, &get_random_bytes(100), n, no_retry);
        let tempvar: f64 = temp.iter().map(|&x| (x as f64 - mu).powi(2)).sum();
        let tempkur: f64 = temp.iter().map(|&x| (x as f64 - mu).powi(4)).sum();

        var += tempvar / (t * 2 * n) as f64;
        kur += tempkur / (t * 2 * n) as f64;
    }

    let end = start.elapsed();

    println!("Sigma for degree {} finished in {:?}", n, end);
    // return data on the following form:
    (mu, var, kur, end)
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
        samples.push(hawksign_x_only(
            &privkey,
            &get_random_bytes(100),
            n,
            no_retry,
        ));
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
