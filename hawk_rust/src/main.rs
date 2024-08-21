use keygen::hawkkeygen;
// use params::initialize_params;
use rngcontext::RngContext;

use crate::utils::{adjoint, bin, int, is_invertible, mod_pow, poly_add, poly_mult_ntt};
use num_bigint::{BigInt, ToBigInt};
use num_traits::{Zero, One};

mod keygen;
mod ntt;
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
    hawkkeygen(8, None);

    // test
    // let a: Vec<BigInt> = bigint_vec(vec![13,4,2,99]);
    // let b: Vec<BigInt> = bigint_vec(vec![8,7, 6, 66]);
    // println!("{:?}", ntru_solve::galois_conjugate(a));
    // let x = (194238123456789098676512312598 as u128).to_bigint().unwrap();
    // let y = (133713371238909804523897987423987 as u128).to_bigint().unwrap();
    // let z = ntru_solve::xgcd(x.clone(), y.clone());
    // println!("gcd({:?}, {:?}) = {:?}",x,y,z);
    //
}

fn print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>())
}

