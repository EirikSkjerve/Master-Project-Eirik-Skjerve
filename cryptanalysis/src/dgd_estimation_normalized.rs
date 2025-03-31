// this file contains experiments for estimating the std. dev. and variance
// for the practical implementation of the "Discrete Gaussian Distribution" in Hawk.

use hawklib::hawkkeygen::hawkkeygen;
use hawklib::hawksign::{hawksign_total, hawksign_x_only};
use hawklib::parameters::{hawk1024_params, hawk256_params, hawk512_params};

use rand::Rng;

use std::fs::File;
use std::io::{stdout, Write};
use std::path::Path;
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};
use std::collections::HashMap;

use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;

use prettytable::{color, Attr, Cell, Row, Table};

pub fn estimate_mem_norm_all(t: usize, store_file: bool) {
    // let ns = vec![256, 512, 1024];
    // let ns = vec![512];
    let ns = vec![256];

    let precision = 20;

    // define table of estimations
    let mut table = Table::new();

    table.add_row(row![
                  i->"Deg",
                  i-> "Num\nVectors",
                  i->"Mu",
                  i->"Var\n(Exp(X^2))",
                  i->"Sigma",
                  i->"Norm.\nvar",
                  i->"Mu4",
                  i->"Mu4-3",
                  i->"Time",
    ]);

    println!("Estimating values for Hawk using {t} samples");
    let start = Instant::now();
    for n in ns {
        // get the right parameters. These are just for comparison
        let (f, sigmasign) = match n {
            256 => (1, hawk256_params::SIGMASIGN),
            512 => (2, hawk512_params::SIGMASIGN),
            _ => (4, hawk1024_params::SIGMASIGN),
        };

        // run the estimation, either sequentially or in parallel
        let (mu, var, normvar, normkur, time) = estimate_mem_norm_par(t, n);
        // let (mu, var, normvar, normkur, time) = estimate_mem_norm(t, n);
        println!("Estimation done for degree {n}");
        // write results to table
        table.add_row(row![
        FG->n.to_string(),
        FW->(t).to_string(),
        FM->format!("{}", mu),
        FB->format!("{}", var),
        Fc->format!("{}", var.sqrt()),
        Fy->format!("{}", normvar),
        FR->format!("{}", normkur),
        Fy->format!("{}", (normkur - 3.0)),
        Fw->format!("{:?}", time)
        ]);
    }
    table.printstd();
    let end = start.elapsed();
    println!("Total time used: {:?}", end);

    // store file
    if store_file {
        let pathname_base = String::from("sigma_table");
        let mut pathname = String::from("sigma_table_0");
        let mut ctr = 1;

        // making sure a file is never written over, and iteratively construct new unique file name
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

pub fn estimate_mem_norm_par(t: usize, n: usize) -> (f64, f64, f64, f64, Duration) {
    let (privkey, _) = hawkkeygen(n);
    let start = Instant::now();

    // let mut min: i64 = 100;
    // let mut max: i64 = -100;
    //
    // for i in 0..t {
    //     let temp: Vec<i64> = hawksign_x_only(&privkey, &get_random_bytes(50), n, true);
    //     for tp in temp {
    //         if tp < min  { min = tp}
    //         if tp > max { max = tp}
    //     }
    // }

    // println!("Min: {min} \nMax: {max}");


    // thread safe variable
    let mut mu: Arc<Mutex<f64>> = Arc::new(Mutex::new(0.0));
    let mut var: Arc<Mutex<f64>> = Arc::new(Mutex::new(0.0));

    let mut frequencies: Arc<Mutex<HashMap<i64, f64>>>= Arc::new(Mutex::new(HashMap::new()));

    // make progressbar with style
    let pb = ProgressBar::new(t as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({per_sec})")
            .unwrap()
            .progress_chars("#>-"),
    );



    // sample x-vectors in parallel
    (0..t).into_par_iter().for_each(|i| {
        let temp: Vec<i64> = hawksign_x_only(&privkey, &get_random_bytes(50), n, false);
        temp.clone().into_iter().for_each(|x| {
            let mut freq_map = frequencies.lock().unwrap();
            *freq_map.entry(x).or_insert(0.0) += 1.0;
        });
        let tempmu: f64 = temp.iter().map(|&x| (x as f64)).sum();

        *mu.lock().unwrap() += tempmu;

        // increment progress bar
        pb.inc(1);
    });

        // Normalize frequencies (convert counts to relative frequencies)
    {
        let mut freq_map = frequencies.lock().unwrap();
        for value in freq_map.values_mut() {
            *value /= (2*t*n) as f64;
        }
    }
        // Sort by key before printing
    let mut sorted_frequencies: Vec<(i64, f64)> = {
        let freq_map = frequencies.lock().unwrap();
        let mut vec: Vec<_> = freq_map.iter().map(|(&k, &v)| (k, v)).collect();
        vec.sort_by_key(|&(k, _)| k); // Sort by key (integer value)
        vec
    };

    // Print sorted relative frequencies
    // for (key, value) in sorted_frequencies {
    //     println!("{}: {:.10}", key, value);
    // }

    let mu = *mu.lock().unwrap() / (2*t*n) as f64;

    pb.finish_with_message("Mu estimation completed");

    // remake progressbar
    let pb = ProgressBar::new(t as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({per_sec})")
            .unwrap()
            .progress_chars("#>-"),
    );

    // sample x-vectors in parallel
    (0..t).into_par_iter().for_each(|i| {
        let temp: Vec<i64> = hawksign_x_only(&privkey, &get_random_bytes(50), n, false);
        // possibly subtract mu from x
        let tempvar: f64 = temp.iter().map(|&x| (x as f64 - mu).powi(2)).sum();

        *var.lock().unwrap() += tempvar;

        // increment progress bar
        pb.inc(1);
    });

    pb.finish_with_message("Sigma^2 estimation completed");

    // remake progressbar
    let pb = ProgressBar::new(t as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({per_sec})")
            .unwrap()
            .progress_chars("#>-"),
    );

    let sigma_sq: f64 = *var.lock().unwrap();
    let sigma = (sigma_sq / (2*t*n) as f64).sqrt();


    // thread-safe variables
    let mut normvar: Arc<Mutex<f64>> = Arc::new(Mutex::new(0.0));
    let mut normkur: Arc<Mutex<f64>> = Arc::new(Mutex::new(0.0));

    // sample x-vectors in parallel
    (0..t).into_par_iter().for_each(|_| {
        let temp: Vec<i64> = hawksign_x_only(&privkey, &get_random_bytes(50), n, false);

        let (tempvar, tempkur): (f64, f64) = temp
            .iter()
            .map(|&x| {
                let normalized = (x as f64 - mu) / sigma;
                // possibly subtract mu from x
                let squared = normalized.powi(2);
                let fourth = normalized.powi(4);
                (squared, fourth)
            })
            .fold((0.0, 0.0), |(sum_var, sum_kur), (squared, fourth)| {
                (sum_var + squared, sum_kur + fourth)
            });

        // update terms
        *normvar.lock().unwrap() += tempvar;
        *normkur.lock().unwrap() += tempkur;

        // increment progress bar
        pb.inc(1);
    });

    pb.finish_with_message("Kurtosis estimation completed");

    let end = start.elapsed();

    // unpack the thread-safe variables into normal variables to return them
    let res_var = *var.lock().unwrap() / (2*t*n) as f64;
    let res_normvar = *normvar.lock().unwrap() / (2*t*n) as f64;
    let res_normkur = *normkur.lock().unwrap() / (2*t*n) as f64;
    (mu, res_var, res_normvar, res_normkur, end)
}

pub fn estimate_mem_norm(t: usize, n: usize) -> (f64, f64, f64, f64, Duration) {
    // create t x-vectors with hawk degree n
    // this function will not store any intermediate samples and therefore will not need a lot of
    // memory. Drawback is that t x-vectors have to be generated thrice; once for mu, once for
    // var and once for normalized kurtosis
    //
    // returns (Exp[x], Var[X], normalized Var[X], normalized Kur[x], time used)
    //

    // generate a keypair
    let (privkey, _) = hawkkeygen(n);

    assert!(n == 256 || n == 512 || n == 1024);

    // for nice printouts
    let mut stdout = stdout();

    // bool determining if sampling of x should retry if a sampled x is not "valid"
    // this has an effect on the measurement of the practical distribution
    let no_retry = true;

    println!(
        "\nEstimating values by sampling {} vectors of length 2*{n}",
        t
    );

    // measure time for mu
    let start = Instant::now();

    // estimating mu
    println!("Using mu=0");
    let mut mu: f64 = 0.0;

    // estimating sigma
    println!("\nEstimating sigma...  ");
    let mut var: f64 = 0.0;
    let mut seedlen = (t / u8::MAX as usize).max(1);
    let mut seed: Vec<u8> = vec![0; seedlen];

    for i in 0..t {
        seed[i % seedlen] += (i & u8::MAX as usize) as u8;
        seed[i % seedlen] %= (u8::MAX);
        // sample x-vectors of length 2n given random "message"
        let temp: Vec<i64> = hawksign_x_only(&privkey, &seed, n, no_retry);
        let tempvar: f64 = temp.iter().map(|&x| (x as f64).powi(2)).sum();
        var += tempvar / (t * 2 * n) as f64;

        // Calculate and display progress
        if i % (t / 100) == 0 || i == t - 1 {
            let progress = (i as f64 / t as f64) * 100.0;
            print!("\rProgress: {:.0}%", progress);
            std::io::Write::flush(&mut std::io::stdout()).unwrap();
        }
    }

    // since variance is sigma squared
    let sigma = var.sqrt();

    println!("");

    // estimate normalized var and normalized kur
    println!("\nEstimating normalized var and kur");
    let (mut normvar, mut normkur): (f64, f64) = (0.0, 0.0);
    let seedlen = (t / u8::MAX as usize).max(1);
    let mut seed: Vec<u8> = vec![0; seedlen];
    for i in 0..t {
        seed[i % seedlen] += 1;
        seed[i % seedlen] %= (u8::MAX);

        let temp: Vec<i64> = hawksign_x_only(&privkey, &seed, n, no_retry);

        // println!("{:?}", temp);
        let (tempvar, tempkur): (f64, f64) = temp
            .iter()
            .map(|&x| {
                let normalized = (x as f64) / sigma - mu;
                let squared = normalized.powi(2);
                let fourth = normalized.powi(4);
                (squared, fourth)
            })
            .fold((0.0, 0.0), |(sum_var, sum_kur), (squared, fourth)| {
                (sum_var + squared, sum_kur + fourth)
            });

        normvar += tempvar / (t * 2 * n) as f64;
        normkur += tempkur / (t * 2 * n) as f64;

        // Calculate and display progress
        if i % (t / 100) == 0 || i == t - 1 {
            let progress = (i as f64 / t as f64) * 100.0;
            print!("\rProgress: {:.0}%", progress);
            std::io::Write::flush(&mut std::io::stdout()).unwrap();
        }
    }

    let end = start.elapsed();

    println!("\nEstimation for degree {} finished in {:?}", n, end);
    // return data on the following form:
    (mu, var, normvar, normkur, end)
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
