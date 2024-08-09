use keygen::hawkkeygen;
use params::initialize_params;
use rngcontext::RngContext;

use crate::utils::{adjoint, bin, int, is_invertible, mod_pow, poly_add, poly_mult_ntt};

mod keygen;
mod ntru_solve;
mod ntt_fft;
mod params;
mod rngcontext;
mod sign;
mod utils;
mod verify;

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
    // initialize_params(8);
    hawkkeygen(8, None);
}

fn print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>())
}
