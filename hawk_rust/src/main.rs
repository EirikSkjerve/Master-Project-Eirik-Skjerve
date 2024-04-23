use rngcontext::RngContext;
use params::initialize_params;

mod keygen;
mod ntru_solve;
mod params;
mod rngcontext;
mod sign;
mod utils;
mod verify;

fn main() {
    let mut rng = RngContext::new(1336);
    let mut rand_bits: u128;
    for _ in 0..10{
        rand_bits = rng.rnd(16);
        println!("Random 16 bits: {}",rand_bits);
    }

    let params = initialize_params(8);
    print_type_of(&params["n"]);
}

fn print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>())
}