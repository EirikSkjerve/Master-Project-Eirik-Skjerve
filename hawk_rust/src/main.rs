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
    // hawkkeygen(8, None);

    let a: Vec<BigInt> = bigint_vec(vec![13,4,2,99]);
    let b: Vec<BigInt> = bigint_vec(vec![8,7, 6, 66]);
    // let c = ntru_solve::karamul(a, b);
    let field_norm_a = ntru_solve::field_norm(&a);
    // println!("{:?} * {:?} mod 2^{} = {:?}", a, b, a.len(), c);
    println!("field norm of {:?} = {:?}", a, field_norm_a);
}

fn print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>())
}

fn bigint_vec(v: Vec<i64>) -> Vec<BigInt> {
    let mut v_big: Vec<BigInt> = Vec::new();
    for i in v.iter() {
        v_big.push(i.to_bigint().unwrap());
    }

    return v_big;
}
