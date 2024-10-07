use keygen::{hawkkeygen_256, hawkkeygen_512, hawkkeygen_1024};
use rngcontext::get_random_bytes;
use sign::{hawksign};
use verify::hawkverify;

use std::time::{Duration, Instant};

mod cdt;
mod codec;
mod compress;
mod decompress;
mod fft;
mod fft_constants;
mod grutils;
mod keygen;
mod ntru_solve;
mod ntt;
mod params;
mod rngcontext;
mod sign;
mod utils;
mod verify;
mod verifyutils;
mod tests;
mod parameters;


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

    // measure time
    let start = Instant::now();
    for _ in 0..samples {

        // generate some random message
        let message = get_random_bytes(100);

        // produce a signature for message
        
        let sig_time = Instant::now();
        let signature = hawksign(logn, &privkey, &message);
        println!("{:?}", sig_time.elapsed());
        
        // let (_salt, _s1) = dec_sig(logn, &signature);

        // verify or reject a message
        let verify = hawkverify(logn, &message, &pubkey, &signature);

        // count failures
        if !verify{
            failed += 1;
        }
    }

    let duration = start.elapsed();
    println!("time used for {} sign and verify on random messages: {:?}", samples, duration);
    println!("{} failed.", failed);
}
