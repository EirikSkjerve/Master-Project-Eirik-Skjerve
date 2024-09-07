use keygen::hawkkeygen;
use sign::sign;
use verify::verify;

// simple rng library
use rand::prelude::*;

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

    loop {
    let mut rng = thread_rng();
    let rand_seed = rng.gen_range(0..99999999);
    // we're generally interested in the lowest security level
    let startkg = Instant::now();
    let keypair = hawkkeygen(8, rand_seed);
    let durkg = startkg.elapsed();

    // public key, secret key
    let (pk, sk) = keypair;
    // println!("pk: {:?} \nsk: {:?}", pk, sk);

    let message = 123456789 as usize;
    // private polynomials in here
    let startsg = Instant::now();
    let signature = sign(8, &sk, message);
    let dursg = startsg.elapsed();
    //
    //

    // public polynomials in here
    let verify = verify(8, message, &pk, &signature);
    if verify {
        break;
    }
    // break;

    }

}

fn test_pipeline() {
    let mut min = 0;
    let mut min_seed = 0;
    let mut sum = 0;

    let mut rng = thread_rng();
    //
    // let randval = rng.gen_range(0..100000000);
    // let key_pair = hawkkeygen(8, randval);
    // let (f, g, F, G, q00, q01, kgseed, counter) = key_pair;
    // let start = Instant::now();
    // for i in 0..10000 {
    //     let message = rng.gen_range(0..10000);
    //     let F_clone = F.clone();
    //     let G_clone = G.clone();
    //     let sig = sign(8, F_clone, G_clone, kgseed, message);
    //     let ver = verify(8, message, q00.clone(), q01.clone(), sig);
    //     if ver {
    //         println!("{}", randval);
    //     }
    // }
    // let duration = start.elapsed();
    // println!("Generated in {:?}", duration);
    // let peak_mem = PEAK_ALLOC.peak_usage_as_kb();
    // println!("Max memory use: {} kb", peak_mem);
}

fn print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>())
}
