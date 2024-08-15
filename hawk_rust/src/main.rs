use keygen::hawkkeygen;
// use params::initialize_params;
use rngcontext::RngContext;

use crate::utils::{adjoint, bin, int, is_invertible, mod_pow, poly_add, poly_mult_ntt};

mod keygen;
mod ntru_solve;
mod ntt_fft;
// mod params;
mod fft;
mod fft_constants;
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

    /*
    // let f: Vec<i64> = vec![1, 1, 0, 0, 0, 0, -1, 2, 0, 0, -1, 1, -2, 1, 0, 0, -2, 1, 1, 1, 0, -1, 0, 0, 0, 0, -1, 1, 1, -2, -1, 2, 1, 0, 0, 1, -1, 1, 1, 1, 0, 1, 1, -1, 0, -1, 1, 0, 2, 1, -1, -1, 0, -1, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, -1, -1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, -1, -1, -1, -2, 0, 0, -2, 0, 1, 0, 0, 1, 1, -2, 0, -1, -1, 0, -1, 0, 1, 1, -1, 0, -1, 0, 0, -1, -1, -1, -1, 0, -1, 1, -2, 1, 1, 2, 1, 1, 0, -1, 0, -2, 0, 0, 0, -1, -2, 1, 0, -1, 1, -1, 0, 0, 0, 1, 1, 2, 1, 0, 1, 0, -1, -2, 0, 0, -1, -1, 1, 1, 1, 1, 0, 0, -1, 1, 1, 1, -1, 2, 1, -1, 0, 1, 1, -1, 1, 0, 0, -1, 0, -2, 1, 0, 1, 1, -1, -1, 0, 1, 0, -1, 0, 0, 0, -2, 0, 2, 0, 0, -2, -1, -2, -1, 1, -1, 1, 1, -1, 2, 0, 0, 2, 0, 0, 0, 0, 0, 0, -1, 1, 2, 0, -1, -1, -1, 1, -1, -1, 1, 0, -1, -1, 0, 0, 0, 0, 1, -1, 0, 2, 2, 1, -1, 1, 0, 1, 0, 0, -1, 0, 0, 1, 0, -2, -2, 1, 1, 0, -2, 0, 2, -1, 1, -1, 1];
    let f: Vec<f64> = vec![1.0, 2.0, 3.0, 4.0];
    let f_fft = fft::fft(&f);
    println!("fft(f): {:?}", f_fft);

    let f_orig = fft::ifft(&f_fft);
    println!("ifft(fft(f)): {:?}", f_orig);

    let p: Vec<i64> = vec![1,2,3,4];
    let inverse_f_fft = fft::inverse_fft(&p);
    println!("inverse f_fft: {:?}", inverse_f_fft);
    */
}

fn print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>())
}
