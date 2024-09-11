use keygen::hawkkeygen;
use rngcontext::get_random_bytes;
use sign::hawksign;
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
    let attempts = 100;
    for i in 0..attempts {
        let logn = 8;

        // let init_seed: [u8;10] = [1,2,3,4,5,6,7,8,9,10];
        let init_seed = get_random_bytes(10);
        let message = get_random_bytes(100);

        let keypair = hawkkeygen(logn, &init_seed);
        let (privkey, pubkey) = keypair;

        let signature = hawksign(logn, &privkey, &message);

        let verify = hawkverify(logn, &message, &pubkey, &signature);
    }

    let duration = start.elapsed();
    println!("time used for 100 keygen, sign and verify: {:?}", duration);
}

fn print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>())
}
