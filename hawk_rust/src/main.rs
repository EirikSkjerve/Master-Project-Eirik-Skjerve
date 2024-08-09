use keygen::hawkkeygen;
use params::initialize_params;
use rngcontext::RngContext;

use crate::utils::{adjoint, bin, int, mod_pow, poly_add, poly_mult};

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

    let h = vec![5, 6, 7, 8];
    // let p = 2147473409;
    let p = 7681;

    let h_hat = ntt_fft::ntt(h, p);
    println!("h hat: {:?}", h_hat);
    let h_orig = ntt_fft::intt(h_hat, p);
    println!("h original: {:?}", h_orig);
}

fn print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>())
}
