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

pub fn estimate_mem_all(t: usize, store_file: bool) {
    // let ns = vec![256, 512, 1024];
    let ns = vec![256];

    let precision = 8;

    // define table of estimations
    let mut table = Table::new();

    table.add_row(row![
                  i->"Deg",
                  i-> "Num\nVectors",
                  i->"Mu",
                  i->"Var (Exp(X^2))",
                  i->"Sigma (Std. Dev)",
                  i->"Th. \nSigma",
                  i->"Diff.",
                  i->"Kur (Exp(X^4))",
                  i->"3 Var ^2",
                  i->"Diff.",
                  i->"Th. \nKur",
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
        let (mu, var, kur, time) = estimate_mem(t, n);

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
        FM->format!("{:.1$}", (3.0*((2.0*sigmasign).powi(4))).abs(), precision),
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

pub fn estimate_mem(t: usize, n: usize) -> (f64, f64, f64, Duration) {
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

    println!("");

    // estimating mu
    let mut mu: f64 = 0.0;
    for i in 0..t {
        stdout.flush().unwrap();
        print!("\r{} %", (((i + 1) as f64 / t as f64) * 100.0).abs() as u8);

        // sample x-vector of length 2n given a random vector
        let temp: i64 = hawksign_x_only(&privkey, &get_random_bytes(100), n, no_retry)
            .iter()
            .sum();
        // mu += temp as f64 / (t * 2 * n) as f64;
        mu += temp as f64
    }

    mu /= (t*2*n) as f64;

    println!("");

    // estimating var and kur
    let mut var: f64 = 0.0;
    let mut kur: f64 = 0.0;

    for i in 0..t {
        stdout.flush().unwrap();
        print!("\r{} %", (((i + 1) as f64 / t as f64) * 100.0).abs() as u8);

        // sample x-vector of length 2n given a random vector
        let temp: Vec<i64> = hawksign_x_only(&privkey, &get_random_bytes(100), n, no_retry);
        let tempvar: f64 = temp.iter().map(|&x| (x as f64 - mu).powi(2)).sum();
        let tempkur: f64 = temp.iter().map(|&x| (x as f64 - mu).powi(4)).sum();

        // var += tempvar / (t * 2 * n) as f64;
        // kur += tempkur / (t * 2 * n) as f64;
        var += tempvar;
        kur += tempkur;
    }

    var /= (t*2*n) as f64;
    kur /= (t*2*n) as f64;

    let end = start.elapsed();

    println!("\nSigma for degree {} finished in {:?}", n, end);
    // return data on the following form:
    (mu, var, kur, end)
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
