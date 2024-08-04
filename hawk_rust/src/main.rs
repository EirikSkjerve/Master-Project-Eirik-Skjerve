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

    println!("{:?}", mod_pow(1337_i32, 999_i32, 150000_i32));
    ntt_fft::primitive_root(54521);
    ntt_fft::primitive_root(2147473409);

}

fn print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>())
}
