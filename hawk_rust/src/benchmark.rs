use hawklib::hawkkeygen::hawkkeygen;
use hawklib::hawksign::hawksign;
use hawklib::hawkverify::hawkverify;

use std::time::Instant;
use rand::Rng;

use prettytable::{Cell, Row, Table};

pub const NUM_SAMPLES: usize = 1000;

pub fn test_all() {
    
    println!("Running tests...");

    // define table of timings
    let mut table = Table::new();
    table.add_row(Row::new(vec![
        Cell::new(&"deg"),
        Cell::new(&"kg (ms)"),
        Cell::new(&"sg (µs)"),
        Cell::new(&"vf (µs)"),
    ]));

    hawkrun(&mut table, 256);
    hawkrun(&mut table, 512);
    hawkrun(&mut table, 1024);

    println!(
        "Average of {} signature generation and verifications",
        NUM_SAMPLES
    );
    table.printstd();
}

pub fn hawkrun(table: &mut Table, n: usize) {
    // start time for key generation
    let kgen_time_start = Instant::now();

    // generate keypair
    let (privkey, pubkey) = hawkkeygen(n);

    // end time for key generation
    let kgen_time_end = kgen_time_start.elapsed();

    // pre-generate some random, uniformly generated messages
    let mut messages: Vec<Vec<u8>> = Vec::with_capacity(NUM_SAMPLES);
    for _ in 0..NUM_SAMPLES {
        messages.push(get_random_bytes(10));
    }

    // create collection of signatures corresponding to messages
    let mut signatures: Vec<(Vec<i64>, Vec<u8>)> = Vec::with_capacity(NUM_SAMPLES);
    // start time for signature generation
    let sig_time_start = Instant::now();
    for i in 0..NUM_SAMPLES {
        // sign the messages
        signatures.push(hawksign(&privkey, &messages[i], n));
    }
    // end time for signature generation
    let sig_time_stop = sig_time_start.elapsed();

    // verify the message/signature pairs
    // start time for signature verification
    let mut num_fails = 0;
    let ver_time_start = Instant::now();
    for i in 0..NUM_SAMPLES {
        let ver = hawkverify(&messages[i], &pubkey, &signatures[i].0, &signatures[i].1, n);
        if !ver {
            num_fails+=1;
        }
    }

    if num_fails != 0 { println!("{} failure(s)", num_fails)}
    // println!("{} failures", num_fails);

    // end time for signature verification
    let ver_time_stop = ver_time_start.elapsed();

    // write all times to table
    table.add_row(Row::new(vec![
        Cell::new(&n.to_string()),
        Cell::new(&kgen_time_end.as_millis().to_string()),
        Cell::new(&(sig_time_stop / NUM_SAMPLES as u32).as_micros().to_string()),
        Cell::new(&(ver_time_stop / NUM_SAMPLES as u32).as_micros().to_string()),
    ]));
}

pub fn get_random_bytes(num_bytes: usize) -> Vec<u8> {
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
