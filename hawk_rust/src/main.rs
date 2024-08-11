use keygen::hawkkeygen;
// use params::initialize_params;
use rngcontext::RngContext;

use crate::utils::{adjoint, bin, int, is_invertible, mod_pow, poly_add, poly_mult_ntt};

mod keygen;
mod ntru_solve;
mod ntt_fft;
// mod params;
mod rngcontext;
mod sign;
mod utils;
mod verify;
mod fft;
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

    // let f: Vec<i64> = vec![1, 1, 0, 0, 0, 0, -1, 2, 0, 0, -1, 1, -2, 1, 0, 0, -2, 1, 1, 1, 0, -1, 0, 0, 0, 0, -1, 1, 1, -2, -1, 2, 1, 0, 0, 1, -1, 1, 1, 1, 0, 1, 1, -1, 0, -1, 1, 0, 2, 1, -1, -1, 0, -1, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, -1, -1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, -1, -1, -1, -2, 0, 0, -2, 0, 1, 0, 0, 1, 1, -2, 0, -1, -1, 0, -1, 0, 1, 1, -1, 0, -1, 0, 0, -1, -1, -1, -1, 0, -1, 1, -2, 1, 1, 2, 1, 1, 0, -1, 0, -2, 0, 0, 0, -1, -2, 1, 0, -1, 1, -1, 0, 0, 0, 1, 1, 2, 1, 0, 1, 0, -1, -2, 0, 0, -1, -1, 1, 1, 1, 1, 0, 0, -1, 1, 1, 1, -1, 2, 1, -1, 0, 1, 1, -1, 1, 0, 0, -1, 0, -2, 1, 0, 1, 1, -1, -1, 0, 1, 0, -1, 0, 0, 0, -2, 0, 2, 0, 0, -2, -1, -2, -1, 1, -1, 1, 1, -1, 2, 0, 0, 2, 0, 0, 0, 0, 0, 0, -1, 1, 2, 0, -1, -1, -1, 1, -1, -1, 1, 0, -1, -1, 0, 0, 0, 0, 1, -1, 0, 2, 2, 1, -1, 1, 0, 1, 0, 0, -1, 0, 0, 1, 0, -2, -2, 1, 1, 0, -2, 0, 2, -1, 1, -1, 1];
    let f: Vec<i64> = vec![1,2,3,4];
    let fft_f = fft::f_fft(&f);
    println!("{:?}", fft_f);

    
    
}

fn print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>())
}
