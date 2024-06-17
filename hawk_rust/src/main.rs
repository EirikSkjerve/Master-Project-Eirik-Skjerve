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

fn main() {
    //hawkkeygen(8, None);
    let a = vec![1, 2];
    let b = vec![4, 3, 2];
    let c = poly_mult(&a, &b);
    println!("{:?} * {:?} = {:?}", a, b, c);
}

fn print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>())
}
