use crate::hawk1024::{
    hawkkeygen_1024::hawkkeygen_1024, hawksign_1024::hawksign_1024,
    hawkverify_1024::hawkverify_1024,
};
use crate::hawk256::{
    hawkkeygen_256::hawkkeygen_256, hawksign_256::hawksign_256, hawkverify_256::hawkverify_256,
};
use crate::hawk512::{
    hawkkeygen_512::hawkkeygen_512, hawksign_512::hawksign_512, hawkverify_512::hawkverify_512,
};

use crate::rngcontext::get_random_bytes;

// need this to reset the table when computing ntt's for different values of n
use crate::ntt_constants::res_z;

use std::time::{Duration, Instant};

use prettytable::{Cell, Row, Table};


pub const NUM_SAMPLES: usize = 500;

pub fn test_all() {
    let mut table = Table::new();
    table.add_row(Row::new(vec![
        Cell::new(&"deg"),
        Cell::new(&"kg (ms)"),
        Cell::new(&"sg (µs)"),
        Cell::new(&"vf (µs)"),
    ]));
    // hawk_256(&mut table);
    // hawk_512(&mut table);
    hawk_1024(&mut table);

    println!("Average of {} signature generation and verifications", NUM_SAMPLES);
    table.printstd();
}

pub fn hawk_256(table: &mut Table) {

    res_z();
    let init_seed = get_random_bytes(15);
    // generate keypair
    // println!("Generating keypair for Hawk {}...", "256".bright_blue());
    let kgen_time_start = Instant::now();
    let (privkey, pubkey) = hawkkeygen_256(&init_seed);
    let kgen_time_end = kgen_time_start.elapsed();
    // println!(
    //     "Hawk 256 keypair generated in {} \n ",
    //     format_duration(kgen_time_end)
    // );
    // println!("Hawk 256 keypair generated in {:?} ", kgen_time_end);

    // pre-generate some messages
    // println!("Generating {} random messages", NUM_SAMPLES);
    let mgen_time_start = Instant::now();
    let mut messages: Vec<Vec<u8>> = Vec::with_capacity(NUM_SAMPLES);
    for _ in 0..NUM_SAMPLES {
        messages.push(get_random_bytes(100));
    }
    let mgen_time_end = mgen_time_start.elapsed();
    // println!(
    //     "Time used generating {} random messages: {:?} \n",
    //     NUM_SAMPLES, mgen_time_end
    // );

    // keep track of number of failed signatures
    let mut num_failed = 0;

    // create collection of signatures corresponding to messages
    let mut signatures: Vec<Vec<u8>> = Vec::with_capacity(NUM_SAMPLES);
    // println!("Generating {} signatures...", NUM_SAMPLES);
    let sig_time_start = Instant::now();
    for i in 0..NUM_SAMPLES {
        signatures.push(hawksign_256(&privkey, &messages[i]));
    }
    let sig_time_stop = sig_time_start.elapsed();
    // println!(
    //     "Time used generating {} signatures: {:?} \n",
    //     NUM_SAMPLES, sig_time_stop
    // );

    // verify the message/signature pairs
    // println!("Verifying {} signatures...", NUM_SAMPLES);
    let ver_time_start = Instant::now();
    for i in 0..NUM_SAMPLES {
        let verification = hawkverify_256(&messages[i], &pubkey, &signatures[i]);

        if !verification {
            num_failed += 1;
        }
    }
    let ver_time_stop = ver_time_start.elapsed();
    // println!(
    //     "Time used verifying {} signatures: {:?}",
    //     NUM_SAMPLES, ver_time_stop
    // );
    table.add_row(Row::new(vec![
        Cell::new("256"),
        Cell::new(&kgen_time_end.as_millis().to_string()),
        Cell::new(&(sig_time_stop / NUM_SAMPLES as u32).as_micros().to_string()),
        Cell::new(&(ver_time_stop / NUM_SAMPLES as u32).as_micros().to_string()),
    ]));
}

pub fn hawk_512(table: &mut Table) {
    res_z();
    let init_seed = get_random_bytes(15);
    // generate keypair
    // println!("Generating keypair for Hawk {}...", "512".bright_blue());
    let kgen_time_start = Instant::now();
    let (privkey, pubkey) = hawkkeygen_512(&init_seed);
    let kgen_time_end = kgen_time_start.elapsed();
    // println!(
    //     "Hawk 512 keypair generated in {} \n ",
    //     format_duration(kgen_time_end)
    // );
    // println!("Hawk 512 keypair generated in {:?} ", kgen_time_end);

    // pre-generate some messages
    // println!("Generating {} random messages", NUM_SAMPLES);
    let mgen_time_start = Instant::now();
    let mut messages: Vec<Vec<u8>> = Vec::with_capacity(NUM_SAMPLES);
    for _ in 0..NUM_SAMPLES {
        messages.push(get_random_bytes(100));
    }
    let mgen_time_end = mgen_time_start.elapsed();
    // println!(
    //     "Time used generating {} random messages: {:?} \n",
    //     NUM_SAMPLES, mgen_time_end
    // );

    // keep track of number of failed signatures
    let mut num_failed = 0;

    // create collection of signatures corresponding to messages
    let mut signatures: Vec<Vec<u8>> = Vec::with_capacity(NUM_SAMPLES);
    // println!("Generating {} signatures...", NUM_SAMPLES);
    let sig_time_start = Instant::now();
    for i in 0..NUM_SAMPLES {
        signatures.push(hawksign_512(&privkey, &messages[i]));
    }
    let sig_time_stop = sig_time_start.elapsed();
    // println!(
    //     "Time used generating {} signatures: {:?} \n",
    //     NUM_SAMPLES, sig_time_stop
    // );

    // verify the message/signature pairs
    // println!("Verifying {} signatures...", NUM_SAMPLES);
    let ver_time_start = Instant::now();
    for i in 0..NUM_SAMPLES {
        let verification = hawkverify_512(&messages[i], &pubkey, &signatures[i]);

        if !verification {
            num_failed += 1;
        }
    }
    let ver_time_stop = ver_time_start.elapsed();
    // println!(
    //     "Time used verifying {} signatures: {:?}",
    //     NUM_SAMPLES, ver_time_stop
    // );
    table.add_row(Row::new(vec![
        Cell::new("512"),
        Cell::new(&kgen_time_end.as_millis().to_string()),
        Cell::new(&(sig_time_stop / NUM_SAMPLES as u32).as_micros().to_string()),
        Cell::new(&(ver_time_stop / NUM_SAMPLES as u32).as_micros().to_string()),
    ]));
}

pub fn hawk_1024(table: &mut Table) {

    res_z();

    let init_seed = get_random_bytes(15);
    // generate keypair
    // println!("Generating keypair for Hawk {}...", "1024".bright_blue());
    let kgen_time_start = Instant::now();
    let (privkey, pubkey) = hawkkeygen_1024(&init_seed);
    let kgen_time_end = kgen_time_start.elapsed();
    // println!(
    //     "Hawk 1024 keypair generated in {} \n ",
    //     format_duration(kgen_time_end)
    // );
    // println!("Hawk 1024 keypair generated in {:?} ", kgen_time_end);

    // pre-generate some messages
    // println!("Generating {} random messages", NUM_SAMPLES);
    let mgen_time_start = Instant::now();
    let mut messages: Vec<Vec<u8>> = Vec::with_capacity(NUM_SAMPLES);
    for _ in 0..NUM_SAMPLES {
        messages.push(get_random_bytes(100));
    }
    let mgen_time_end = mgen_time_start.elapsed();
    // println!(
    //     "Time used generating {} random messages: {:?} \n",
    //     NUM_SAMPLES, mgen_time_end
    // );

    // keep track of number of failed signatures
    let mut num_failed = 0;

    // create collection of signatures corresponding to messages
    let mut signatures: Vec<Vec<u8>> = Vec::with_capacity(NUM_SAMPLES);
    // println!("Generating {} signatures...", NUM_SAMPLES);
    let sig_time_start = Instant::now();
    for i in 0..NUM_SAMPLES {
        signatures.push(hawksign_1024(&privkey, &messages[i]));
    }
    let sig_time_stop = sig_time_start.elapsed();
    // println!(
    //     "Time used generating {} signatures: {:?} \n",
    //     NUM_SAMPLES, sig_time_stop
    // );
    //
    // // verify the message/signature pairs
    // println!("Verifying {} signatures...", NUM_SAMPLES);
    let ver_time_start = Instant::now();
    for i in 0..NUM_SAMPLES {
        let verification = hawkverify_1024(&messages[i], &pubkey, &signatures[i]);

        if !verification {
            num_failed += 1;
        }
    }
    let ver_time_stop = ver_time_start.elapsed();
    // println!(
    //     "Time used verifying {} signatures: {:?}",
    //     NUM_SAMPLES, ver_time_stop
    // );
    table.add_row(Row::new(vec![
        Cell::new("1024"),
        Cell::new(&kgen_time_end.as_millis().to_string()),
        Cell::new(&(sig_time_stop / NUM_SAMPLES as u32).as_micros().to_string()),
        Cell::new(&(ver_time_stop / NUM_SAMPLES as u32).as_micros().to_string()),
    ]));

}
