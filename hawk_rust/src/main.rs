use keygen::hawkkeygen;
use rngcontext::RngContext;
use sign::sign;
use verify::verify;

// simple rng library
use rand::prelude::*;

use crate::utils::{
    adjoint, bin, int, is_invertible, mod_pow, modulo, poly_add, poly_mult_ntt, poly_sub,
};
use num_bigint::{BigInt, ToBigInt};
use num_traits::{One, Zero};

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
    let mut rng = thread_rng();
    let rand_seed = rng.gen_range(0..99999999);
    // we're generally interested in the lowest security level
    let startkg = Instant::now();
    let keypair = hawkkeygen(8, 123222);
    let durkg = startkg.elapsed();

    test_compress();
    // let message = 123456789 as usize;
    // // private polynomials in here
    // let startsg = Instant::now();
    // let signature = sign(8, F, G, kgseed, message);
    // let dursg = startsg.elapsed();
    //
    // println!("Keygen: {:?}", durkg);
    // println!("Signature: {:?}", dursg);
    //
    // // public polynomials in here
    // let verify = verify(8, message, q00, q01, signature);
    // println!("verify: {}", verify);
    //
}

fn test_compress() {
    let a: Vec<i64> = vec![0;256];
    let b: Vec<i64> = vec![1;256];
    let enc_pub = codec::enc_pub(8, &a, &b);
    println!("enc pub: {:?}", enc_pub);
    let ab = codec::dec_pub(8, &enc_pub);
    println!("original: {:?}", ab);
}

fn test_func() {
    let a: Vec<i32> = vec![1, 2, 3, 4, 5, 6, 7, 8];
    let b: Vec<i32> = vec![-2, -4, -6, -1];
    let a_ifft = verifyutils::ifft(&a);
    println!("a invfft: {:?}", a_ifft);

    // println!("a_fft = {:?} \nb_fft = {:?}", a_fft, b_fft);
}

fn test_pipeline() {
    let mut min = 0;
    let mut min_seed = 0;
    let mut sum = 0;

    let mut rng = thread_rng();

    let randval = rng.gen_range(0..100000000);
    let key_pair = hawkkeygen(8, randval);
    let (f, g, F, G, q00, q01, kgseed, counter) = key_pair;
    let start = Instant::now();
    for i in 0..10000 {
        let message = rng.gen_range(0..10000);
        let F_clone = F.clone();
        let G_clone = G.clone();
        let sig = sign(8, F_clone, G_clone, kgseed, message);
        let ver = verify(8, message, q00.clone(), q01.clone(), sig);
        if ver {
            println!("{}", randval);
        }
    }
    let duration = start.elapsed();
    println!("Generated in {:?}", duration);
    let peak_mem = PEAK_ALLOC.peak_usage_as_kb();
    println!("Max memory use: {} kb", peak_mem);
}

fn print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>())
}
