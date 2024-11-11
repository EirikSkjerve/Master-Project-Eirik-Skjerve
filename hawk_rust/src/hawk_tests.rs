use crate::hawkkeygen::hawkkeygen;
use crate::hawksign::hawksign;
use crate::hawkverify::hawkverify;

use crate::rngcontext::get_random_bytes;

use std::time::{Duration, Instant};

use prettytable::{Cell, Row, Table};

pub const NUM_SAMPLES: usize = 1000;

pub fn test_all() {
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
    // generate keypair
    let kgen_time_start = Instant::now();
    let (privkey, pubkey) = hawkkeygen(n);
    let kgen_time_end = kgen_time_start.elapsed();

    // pre-generate some messages
    let mut messages: Vec<Vec<u8>> = Vec::with_capacity(NUM_SAMPLES);
    for _ in 0..NUM_SAMPLES {
        messages.push(get_random_bytes(100));
    }

    // keep track of number of failed signatures
    let mut num_failed = 0;

    // create collection of signatures corresponding to messages
    let mut signatures: Vec<(Vec<i64>, Vec<u8>)> = Vec::with_capacity(NUM_SAMPLES);
    let sig_time_start = Instant::now();
    for i in 0..NUM_SAMPLES {
        signatures.push(hawksign(&privkey, &messages[i], n));
    }
    let sig_time_stop = sig_time_start.elapsed();

    // verify the message/signature pairs
    let ver_time_start = Instant::now();
    for i in 0..NUM_SAMPLES {
        let verification = hawkverify(&messages[i], &pubkey, &signatures[i].0, &signatures[i].1, n);

        if !verification {
            num_failed += 1;
        }
    }

    let ver_time_stop = ver_time_start.elapsed();
    table.add_row(Row::new(vec![
        Cell::new(&n.to_string()),
        Cell::new(&kgen_time_end.as_millis().to_string()),
        Cell::new(&(sig_time_stop / NUM_SAMPLES as u32).as_micros().to_string()),
        Cell::new(&(ver_time_stop / NUM_SAMPLES as u32).as_micros().to_string()),
    ]));
}
