use keygen::hawkkeygen;
use params::initialize_params;
use rngcontext::RngContext;

use crate::utils::{adjoint, bin, int, mod_pow, poly_add, poly_mult_ntt, is_invertible};

mod keygen;
mod ntru_solve;
mod ntt_fft;
mod params;
mod rngcontext;
mod sign;
mod utils;
mod verify;

fn main() {
    // initialize_params(8);
    // hawkkeygen(8, None);

    /*
    let test = ntt_fft::get_roots(2147473409, 256);
    let zetas = test.0;
    let izetas = test.1;
    println!("zetas: {:?}",zetas);
    println!("izetas: {:?}",izetas);
    */

    let f = vec![5, 6, 7, 8];
    let g = vec![1, 2, 3, 4];

    // let p = 2147473409;
    let p = 7681;

    //let fg = poly_mult_ntt(f, g, p);
    //println!("f*g = {:?}", fg);
    let invertible = is_invertible(&f, p);
    println!("{}", invertible);
}

fn print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>())
}
