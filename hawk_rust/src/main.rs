use keygen::hawkkeygen;
// use params::initialize_params;
use rngcontext::RngContext;

use crate::utils::{adjoint, bin, int, is_invertible, mod_pow, poly_add, poly_mult_ntt};

mod keygen;
mod ntt_fft;
// mod params;
mod fft;
mod fft_constants;
mod rngcontext;
mod sign;
mod utils;
mod verify;
mod ntru_solve;
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
    // hawkkeygen(8, None);

    let a = vec![13,4,2,99];
    let b = vec![8,7, 6, 66];
    let c = ntru_solve::karamul(&a, &b);
    println!("{:?} * {:?} mod 2^{} = {:?}", a, b, a.len(), c);
}

fn print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>())
}
