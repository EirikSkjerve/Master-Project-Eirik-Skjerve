use keygen::hawkkeygen;
use params::initialize_params;
use rngcontext::RngContext;

use crate::utils::{bin, int};

mod keygen;
mod ntru_solve;
mod params;
mod rngcontext;
mod sign;
mod utils;
mod verify;

fn main() {
    hawkkeygen(8, None);
}

fn print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>())
}
