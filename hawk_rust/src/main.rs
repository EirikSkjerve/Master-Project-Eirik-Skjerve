
use keygen::hawkkeygen;
use params::initialize_params;
use rngcontext::RngContext;

use crate::utils::{adjoint, bin, int, poly_add, poly_mult, mod_pow};

mod keygen;
mod ntru_solve;
mod params;
mod rngcontext;
mod sign;
mod utils;
mod verify;
mod ntt_fft;

fn main() {
    initialize_params(8);
    hawkkeygen(8, None);

    let test = ntt_fft::get_roots(2147473409, 256);
    let zetas = test.0;
    let izetas = test.1;
    println!("zetas: {:?}",zetas);
}

fn print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>())
}
