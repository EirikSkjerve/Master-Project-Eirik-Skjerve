use rngcontext::RngContext;

mod keygen;
mod ntru_solve;
mod params;
mod rngcontext;
mod sign;
mod utils;
mod verify;

fn main() {
    let mut rng = RngContext::new(1336);
    let mut rand_bits: u128 = 0;
    for i in 0..10{
        rand_bits = rng.rnd(16);
        println!("Random 16 bits: {}",rand_bits);
    }
}
