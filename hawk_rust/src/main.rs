use keygen::hawkkeygen;
use sign::sign;
use rngcontext::RngContext;

// simple rng library
use rand::prelude::*;

use crate::utils::{adjoint, bin, int, is_invertible, mod_pow, poly_add, poly_mult_ntt};
use num_bigint::{BigInt, ToBigInt};
use num_traits::{One, Zero};

use std::time::{Duration, Instant};

mod keygen;
mod ntt;
// mod params;
mod fft;
mod fft_constants;
mod ntru_solve;
mod rngcontext;
mod sign;
mod utils;
mod verify;
mod cdt;

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
    let mut rng = thread_rng();
    let rand_seed = rng.gen_range(0..99999999);

    // we're generally interested in the lowest security level
    let keypair = hawkkeygen(8, 123222);
    let (f, g, F, G, q00, q01, kgseed, counter) = keypair;

    let signature = sign(8, F, G, kgseed, 123456789);
    // test_pipeline();
}

fn test_pipeline() {
    let mut min = 0;
    let mut min_seed = 0;
    let mut sum = 0;

    let mut rng = thread_rng();


    let start = Instant::now();
    for i in 0..100{
        let randval = rng.gen_range(0..100000000);
        println!("iteration: {} random value: {}", i, randval);
        let key_pair = hawkkeygen(8, randval);
        let counter = key_pair.7;
        sum += counter;
        if counter < min{
            min = counter;
            min_seed = i;
        }
    }
    println!("Average attempts: {}",sum as f64/100.0);
    println!("attempts: {} on seed {}", min, min_seed);
    let duration = start.elapsed();
    // println!("Keys: {:?}", key_pair);
    println!("Generated in {:?} s", duration);
    let peak_mem = PEAK_ALLOC.peak_usage_as_kb();
    println!("Max memory use: {} kb", peak_mem);


}

fn print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>())
}
