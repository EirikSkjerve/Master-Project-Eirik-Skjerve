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
        // let keypair = hawkkeygen(8, rand_seed);
        let keypair = hawkkeygen(8, 1234);
        let durkg = startkg.elapsed();

        // public key, secret key
        let (pk, sk) = keypair;
        // println!("pk: {:?} \nsk: {:?}", pk, sk);

        let message = "Test message 2";
        // private polynomials in here
        let startsg = Instant::now();
        let signature = sign(8, &sk, message);
        // println!("sig1: {:?} \nsig2: {:?}", signature, signature_2);
        // println!("signature: {:?}", signature);
        let dursg = startsg.elapsed();
        //
        //
        // public polynomials in here
        let verify = verify(8, message, &pk, &signature);
        if verify {
            break;
        }
        break;
    }
}

fn print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>())
}
