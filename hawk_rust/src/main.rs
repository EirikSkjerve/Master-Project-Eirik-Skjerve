use rngcontext::get_random_bytes;

use std::time::{Duration, Instant};

mod fft;
mod fft_constants;
mod grutils;
mod ntru_solve;
mod ntt;
mod params;
mod rngcontext;
mod utils;
mod verifyutils;
mod tests;

mod parameters;

mod hawk256;
mod hawk512;

use hawk256::{hawkkeygen_256::hawkkeygen_256, hawksign_256::hawksign_256, hawkverify_256::hawkverify_256};

// mod hawk1024;

mod compression;

// memory measurement
use peak_alloc::PeakAlloc;

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

/*
   Driver code for HAWK implementation.

   Integer types:
   n: 256-1024, use u16,
   log n: 8-10, use u8,
   vectors of bits: 0/1, use u8,
   vectors f and g: initially low numbers, but is changed fast.
       set as Vec<i64>, will have negative values

*/

fn main() {
    test1();
}

fn test1() {

    // number of samples
    let samples = 1;

    // set HAWK-degree
    let logn = 8;
    // some initial seed for the key generation process
    let init_seed = get_random_bytes(10);
    // the computed keypair
    let keypair = hawkkeygen_256(&init_seed);
    // unpack the keys
    let (privkey, pubkey) = &keypair;

    // keepng track of failed signatures
    let mut failed = 0;

    for _ in 0..samples {

        // generate some random message
        let message = get_random_bytes(100);

        // produce a signature for message
        let sig_time = Instant::now();
        let signature = hawksign_256(&privkey, &message);
        println!("Sig time {:?}", sig_time.elapsed());
        
        // verify or reject a message
        let ver_time = Instant::now();
        let verify = hawkverify_256(&message, &pubkey, &signature);
        println!("Ver time {:?}",ver_time.elapsed());

        // count failures
        if !verify{
            failed += 1;
        }
    }

}
