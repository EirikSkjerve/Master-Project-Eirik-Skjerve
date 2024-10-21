use nalgebra as na;
use na::*;

use crate::rngcontext::get_random_bytes;

const NUM_SAMPLES: usize = 100;
const N: usize = 10;

fn get_uni_slice_int(n: usize, entry_bound: usize) -> Vec<i16> {

    // get some random bytes (i.e. integers in [0..255], that are uniformly distributed
    let rnd_bytes_u8: Vec<u8> = get_random_bytes(n);

    // make copy with reduced entries mod n, and convert them to i16 for possible negation
    let mut rnd_bytes_i16: Vec<i16> = rnd_bytes_u8.iter().map(|&x| (x % (entry_bound+1) as u8) as i16).collect();
    for (i, rb) in rnd_bytes_u8.iter().enumerate(){
        // negating 50% of the entries
        if *rb >= 128 {
            rnd_bytes_i16[i] = -rnd_bytes_i16[i];
        }
    }
    rnd_bytes_i16

}

fn get_uni_slice_float(n: usize, dist_bound: usize) -> Vec<f64> {
    
    // get some random bytes, uniformly distributed
    let rnd_bytes_u8: Vec<u8> = get_random_bytes(n);

    // make a copy with normalised entries
    let mut rnd_bytes_f64: Vec<f64> = rnd_bytes_u8
        .iter()
        .map(|&x| (x as f64) / (255.0 / dist_bound as f64))
        .collect();

    for (i, rb) in rnd_bytes_u8.iter().enumerate(){
        // negating 50% of the entries
        if *rb >= 128 {
            rnd_bytes_f64[i] = -rnd_bytes_f64[i];
        }
    }
    rnd_bytes_f64
}

/// generate a secret matrix V than we will try and retrieve
fn gen_sec_mat(degree: usize, entry_bound: usize) -> Matrix<i16, Dyn, Dyn, VecStorage<i16, Dyn, Dyn>>{

    // generate slice of uniformly random integers between -entry_bound..entry_bound
    let uni_slice = get_uni_slice_int(degree*degree, entry_bound);

    // create a matrix from this slice
    let sec_v = DMatrix::from_row_slice(degree, degree, &uni_slice);
    println!("{:?}", sec_v);
    println!("element at 0, 1: {}", sec_v[(0,1)]);

    // return secret matrix v
    sec_v
}


pub fn run_hpp_attack() {

    let entry_bound = 2;
    let dist_bound = 1;

    // create random secret vector v
    let sec_v = gen_sec_mat(N, entry_bound);
    // create a copy of v with floating point numbers
    let sec_v_f = sec_v.map(|x| x as f64); 
    // just hardcoding this for now
    let ex2 = 0.333;

    let mut uni_samples = vec![];

    for i in 0..NUM_SAMPLES {
        let x = get_uni_slice_float(N, dist_bound);
        let x_vec = DVector::from_row_slice(&x);
        let y_vec = x_vec.transpose()*sec_v_f.clone();
        uni_samples.push(y_vec);
    }

}
