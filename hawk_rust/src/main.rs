use keygen::hawkkeygen;
use params::initialize_params;
use rngcontext::RngContext;

use crate::utils::{bin, int, adjoint};

mod keygen;
mod ntru_solve;
mod params;
mod rngcontext;
mod sign;
mod utils;
mod verify;

fn main() {
    hawkkeygen(8, None);
    let test = vec![-2,1,0,0,1,-1,2,0];
    let adj_test = adjoint(&test);
    println!("{:?}", adj_test);
}

fn print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>())
}
