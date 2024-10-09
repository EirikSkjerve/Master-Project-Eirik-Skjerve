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
use std::time::Instant;

pub const NUM_SAMPLES: usize = 1000;

pub fn test_all() {
    hawk_256();
    // hawk_512();
    // hawk_1024();
}

pub fn hawk_256() {
    let init_seed = get_random_bytes(15);

    // generate keypair
    println!("Generating keypair for Hawk 256...");
    let kgen_time_start = Instant::now();
    let (privkey, pubkey) = hawkkeygen_256(&init_seed);
    let kgen_time_end = kgen_time_start.elapsed();
    println!("Hawk 256 keypair generated in {:?} \n", kgen_time_end);

    // pre-generate some messages
    println!("Generating {} random messages", NUM_SAMPLES);
    let mgen_time_start = Instant::now();
    let mut messages: Vec<Vec<u8>> = Vec::with_capacity(NUM_SAMPLES);
    for _ in 0..NUM_SAMPLES {
        messages.push(get_random_bytes(100));
    }
    let mgen_time_end = mgen_time_start.elapsed();
    println!(
        "Time used generating {} random messages: {:?} \n",
        NUM_SAMPLES, mgen_time_end
    );

    let mut num_failed = 0;

    let mut signatures: Vec<Vec<u8>> = Vec::with_capacity(NUM_SAMPLES);

    println!("Generating {} signatures...", NUM_SAMPLES);
    let sig_time_start = Instant::now();
    for i in 0..NUM_SAMPLES {
        signatures.push(hawksign_256(&privkey, &messages[i]));
    }
    let sig_time_stop = sig_time_start.elapsed();
    println!(
        "Time used generating {} signatures: {:?} \n",
        NUM_SAMPLES, sig_time_stop
    );

    println!("Verifying {} signatures...", NUM_SAMPLES);
    let ver_time_start = Instant::now();
    for i in 0..NUM_SAMPLES {
        let verification = hawkverify_256(&messages[i], &pubkey, &signatures[i]);

        if !verification {
            num_failed += 1;
        }
    }
    let ver_time_stop = ver_time_start.elapsed();
    println!(
        "Time used verifying {} signatures: {:?}",
        NUM_SAMPLES, ver_time_stop
    );
    println!("Number of verifications failed: {} \n", num_failed);
}

pub fn hawk_512() {
    let init_seed = get_random_bytes(15);

    // generate keypair
    let (privkey, pubkey) = hawkkeygen_512(&init_seed);

    let mut num_failed = 0;

    for i in 0..NUM_SAMPLES {
        let message = get_random_bytes(100);

        let sig_time_start = Instant::now();
        let signature = hawksign_512(&privkey, &message);
        let sig_time_stop = sig_time_start.elapsed();

        let ver_time_start = Instant::now();
        let verification = hawkverify_512(&message, &pubkey, &signature);
        let ver_time_stop = ver_time_start.elapsed();

        if !verification {
            num_failed += 1;
        }
    }
}

pub fn hawk_1024() {
    let init_seed = get_random_bytes(15);

    // generate keypair
    let (privkey, pubkey) = hawkkeygen_1024(&init_seed);

    let mut num_failed = 0;

    for i in 0..NUM_SAMPLES {
        let message = get_random_bytes(100);

        let sig_time_start = Instant::now();
        let signature = hawksign_1024(&privkey, &message);
        let sig_time_stop = sig_time_start.elapsed();

        let ver_time_start = Instant::now();
        let verification = hawkverify_1024(&message, &pubkey, &signature);
        let ver_time_stop = ver_time_start.elapsed();

        if !verification {
            num_failed += 1;
        }
    }
}
