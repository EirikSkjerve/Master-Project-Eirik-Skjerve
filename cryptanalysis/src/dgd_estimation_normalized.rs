// this file contains experiments for estimating the std. dev. and variance
// for the practical implementation of the "Discrete Gaussian Distribution" in Hawk.

use hawklib::hawkkeygen::hawkkeygen;
use hawklib::hawksign::{hawksign_total, hawksign_x_only};
use hawklib::parameters::{hawk1024_params, hawk256_params, hawk512_params};

use rand::Rng;
use std::fs::File;
use std::io::{stdout, Write};
use std::path::Path;
use std::time::{Duration, Instant};

use prettytable::{color, Attr, Cell, Row, Table};

pub fn estimate_mem_norm_all(t: usize, store_file: bool) {
    // let ns = vec![256, 512, 1024];
    let ns = vec![512];

    let precision = 8;

    // define table of estimations
    let mut table = Table::new();

    table.add_row(row![
                  i->"Deg",
                  i-> "Num\nVectors",
                  i->"Mu",
                  i->"Var\n(Exp(X^2))",
                  i->"Est. Sigma",
                  i->"Norm.\nvar",
                  i->"Mu4",
                  i->"|3-Mu4|",
                  i->"Time",
    ]);

    for n in ns {
        // get the right parameters. These are just for comparison
        let (f, sigmasign) = match n {
            256 => (1, hawk256_params::SIGMASIGN),
            512 => (2, hawk512_params::SIGMASIGN),
            _ => (4, hawk1024_params::SIGMASIGN),
        };

        // run the estimation
        let (mu, var, normvar, normkur, time) = estimate_mem_norm(t / f, n);

        // write results to table
        table.add_row(row![
        FG->n.to_string(),
        FW->(t/f).to_string(),
        FM->format!("{}", mu),
        FB->format!("{:.1$}", var, precision),
        Fc->format!("{:.1$}", var.sqrt(), precision),
        Fy->format!("{}", normvar),
        FR->format!("{:.1$}", normkur, precision),
        Fy->format!("{:.1$}", (3.0 - normkur).abs(), precision),
        Fw->format!("{:?}", time)
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

pub fn estimate_mem_norm(t: usize, n: usize) -> (f64, f64, f64, f64, Duration) {
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
    let mut seedlen = (t/u8::MAX as usize).max(1);
    let mut seed: Vec<u8> = vec![0;seedlen];
    for i in 0..t {
        seed[i % seedlen] += 1;
        seed[i % seedlen] %= (u8::MAX);
        // sample x-vectors of length 2n given random "message"
        let temp: Vec<i64> = hawksign_x_only(&privkey, &seed, n, no_retry);
        // println!("{:?}", temp);
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
    let seedlen = (t/u8::MAX as usize).max(1);
    let mut seed: Vec<u8> = vec![0;seedlen];
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

        // let tempkur: f64 = temp
        //     .iter()
        //     .map(|&x| (x as f64/sigma).powi(4) )
        //     .sum();
        //
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
