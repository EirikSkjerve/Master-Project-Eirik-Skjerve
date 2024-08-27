use sha3::{
    digest::{ExtendableOutput, ExtendableOutputReset, Update},
    Shake256,
};

use rand::{rngs::StdRng, Rng, SeedableRng};
use rand::prelude::*;


use crate::cdt::get_table;
use crate::keygen::generate_f_g;
use crate::rngcontext::RngContext;
use crate::utils::{bytes_to_poly, modulo, poly_add, poly_mult_ntt};

pub fn sample(seed: Vec<u8>, t: Vec<i64>, n: u16) {
    let (T0, T1) = get_table();

}

pub fn sign(logn: u8, F: Vec<i64>, G: Vec<i64>, kgseed: usize, msg: usize) {

    let mut Fmod2: Vec<i64> = Vec::with_capacity(F.len());
    let mut Gmod2: Vec<i64> = Vec::with_capacity(G.len());

    for i in 0..F.len() {
        Fmod2.push(modulo(F[i],2));
        Gmod2.push(modulo(G[i],2));
    }
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

        let (h0, h1) = (&bytes_to_poly(&h[0..256/8], n), &bytes_to_poly(&h[(256/8)..256/4], n));
        // convert h to two polynomials
        
        // compute target vectors t0, t1
            

        let mut t0: Vec<u8> = Vec::with_capacity(n);
        let mut t1: Vec<u8> = Vec::with_capacity(n);

        let temp_t0 = poly_add(&poly_mult_ntt(&h0, &f, p), &poly_mult_ntt(&h1, &Fmod2, p));
        let temp_t1 = poly_add(&poly_mult_ntt(&h0, &g, p), &poly_mult_ntt(&h1, &Gmod2, p));

        for i in 0..n {
            // we can be sure these values fit inside an u8 since values are 0 and 1
            t0.push(modulo(temp_t0[i], 2) as u8);
            t1.push(modulo(temp_t1[i], 2) as u8);
        }

        // get random seed = M || kgseed || a+1 || rnd(320)
        let seed = rng.rnd(40);

        // compute (x0, x1) from sample()

        // continue if some requirements are not fulfilled 
        
        // compute and return sig = salt, s1
        break;
    }
}

fn add_bytes(arr: &[&[u8]]) -> Vec<u64> {
    let mut n_temp = 0;
    for a in arr.iter(){
        let a_max = a.iter().max().unwrap();
        if *a_max > n_temp {
            n_temp = *a_max;
        }
    }

    let n = n_temp as usize;

    let mut res: Vec<u64> = Vec::with_capacity(n);
    let mut temp = 0;
    for i in 0..n {
        temp = 0;
        for a in arr.iter(){
            if i < a.len(){
                temp += a[i]; 
            }
        }
        res.push(temp as u64);
    }
    return res;
}
