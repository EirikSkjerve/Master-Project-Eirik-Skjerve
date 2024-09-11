use keygen::hawkkeygen;
use rngcontext::get_random_bytes;
use sign::hawksign;
use verify::hawkverify;
use codec::dec_sig;

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
    let start = Instant::now();
    let samples = 1000;

    let logn = 8;
    let init_seed = get_random_bytes(10);
    let keypair = hawkkeygen(logn, &init_seed);
    let (privkey, pubkey) = &keypair;

    let mut failed = 0;
    for _ in 0..samples {

        let message = get_random_bytes(100);

        let signature = hawksign(logn, &privkey, &message);
        let (_salt, _s1) = dec_sig(logn, &signature);

        let verify = hawkverify(logn, &message, &pubkey, &signature);
        if !verify{
            failed += 1;
        }
    }

    let duration = start.elapsed();
    println!("time used for {} sign and verify on random messages: {:?}", samples, duration);
    println!("{} failed.", failed);
}
