use sha3::{
    digest::{ExtendableOutput, ExtendableOutputReset, Update},
    Shake256,
};

use rand::{rngs::StdRng, Rng, SeedableRng};
use rand::prelude::*;


use crate::cdt::get_table;
use crate::keygen::generate_f_g;
use crate::rngcontext::RngContext;

pub fn sample(seed: Vec<u8>, t: Vec<i64>, n: u16) {
    let (T0, T1) = get_table();

}

pub fn sign(logn: u8, F: Vec<i64>, G: Vec<i64>, kgseed: usize, msg: usize) {

    // initialize a new RngContext with some random seed
    // use random() instead of fixed seed
    let mut rng = RngContext::new(1337);

    let n = 1 << logn;
    let (f, g) = generate_f_g(kgseed, logn);

    // compute hash M
    let mut shaker = Shake256::default();
    shaker.update(&msg.to_ne_bytes());
    let mut m: [u8; 64] = [0; 64];
    shaker.finalize_xof_reset_into(&mut m);
    

    let mut a: usize = 0;
    let p = (1<<16) + 1;

    loop{
        // compute salt
        shaker.update(&m);
        shaker.update(&kgseed.to_ne_bytes());
        shaker.update(&a.to_ne_bytes());
        shaker.update(&rng.rnd(14).to_ne_bytes());
        let mut salt: [u8; 14] = [0; 14];
        // resets the hasher instance
        shaker.finalize_xof_reset_into(&mut salt);
        println!("Salt: {:?}", salt);
        
        // compute new hash h

        shaker.update(&m);
        shaker.update(&salt);
        // the 256 is the degree. Should depend on input logn
        let mut h: [u8; 256/4] = [0; 256/4];
        shaker.finalize_xof_reset_into(&mut h);

        // convert h to two polynomials
        
        // compute target vectors t0, t1

        // get random seed = M || kgseed || a+1 || rnd(320)

        // compute (x0, x1) from sample()

        // continue if some requirements are not fulfilled 
        
        // compute and return sig = salt, s1
        break;
    }
}
