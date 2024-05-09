use keygen::hawkkeygen;
use rngcontext::RngContext;
use params::initialize_params;

use crate::utils::{bin, int};

mod keygen;
mod ntru_solve;
mod params;
mod rngcontext;
mod sign;
mod utils;
mod verify;

fn main() {

    // testing rngContext
    let mut rng = RngContext::new(1336);
    let mut rand_bits: u128;
    for _ in 0..10{
        rand_bits = rng.rnd(16);
        println!("Random 16 bits: {}",rand_bits);
    }

    hawkkeygen(8, None);

}

fn print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>())
}