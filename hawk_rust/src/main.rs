use keygen::hawkkeygen;
use params::initialize_params;
use rngcontext::RngContext;

use crate::utils::{adjoint, bin, int, poly_add, poly_mult};

mod keygen;
mod ntru_solve;
mod params;
mod rngcontext;
mod sign;
mod utils;
mod verify;
mod ntt;

fn main() {
    initialize_params(8);
    hawkkeygen(8, None);
}

fn print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>())
}
