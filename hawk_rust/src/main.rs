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

    // testing int / bin functions
    let a: u128 = 31;
    let a_bin:[u8;5] = [1,1,1,1,1];
    let a_bin_vec: Vec<u8> = vec![1,1,1,1,1];

    println!("{:?}", bin(a, 7));
    println!("{}", int(a_bin));
    println!("{}", int(a_bin_vec));

}

fn print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>())
}