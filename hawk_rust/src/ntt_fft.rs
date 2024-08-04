// This is the Number Theoretic Transform module
// see https://eprint.iacr.org/2024/585.pdf

use num::One;
use prime_factorization::Factorization;
use crate::utils::{mod_pow, int, bin};
use num_bigint::{BigUint, ToBigUint};

// ntt of a polynomial
pub fn ntt(f: Vec<i32>, p: u32) -> Vec<i32> {

    let mut ntt_f = vec![0;f.len()];

    // using *this* algorithm

    return ntt_f;
}

// implementation of bit reversal of an integer
pub fn brv(x: u128, log_n: u128) -> u128 {
    // following is HAWK-teams implementation. Might go about it a different way
    // no idea how it works though

    let mut r = 0;
    for i in 0..log_n{
        r |= ((x >> i) & 1) << (log_n - 1 - i);
    }
    return r;
}

// given prime p and order n, compute 2n-th root of unity mod p
pub fn get_roots(p: u128, n: u128) -> (Vec<u128>, Vec<u128>) {


    let mut g0 = primitive_root(p);
    let b = (p-1)/(2*n);
    g0 = mod_pow(g0, b, p);
    println!("{}", g0);

    let zetas = compute_zetas(g0, p, n);    
    let mut izetas: Vec<u128> = Vec::new();

    for z in zetas.iter() {
        izetas.push(mod_pow(*z, p-2, p));
    }

    return (zetas, izetas);
}

// comptute zeta values in a specific order for usage in ntt/intt functions
fn compute_zetas(root: u128, p: u128, n: u128) -> Vec<u128>{

    let mut zetas: Vec<u128> = Vec::new();
    let log_n = (n as f64).log2() as u128;
    println!("logn: {}", log_n);

    for u in 0.. n {

        let brv = brv(u, log_n);
        //println!("brv: {}", brv);
        let res = mod_pow(root, brv, p);
        //println!("res: {}", res);
        zetas.push(mod_pow(root, brv, p));
    }
    println!("{}", zetas.len());
    return zetas;
}

// return the primitive root of prime p
// using https://www.geeksforgeeks.org/primitive-root-of-a-prime-number-n-modulo-n/
pub fn primitive_root(p: u128) -> u128 {

    let mut g: u128 = 2;

    // phi(p) = p-1
    let s = p - 1;

    // compute prime factors of p-1
    let mut s_factors = Factorization::run(s).factors;

    // remove duplicates
    s_factors.dedup();

    // check if g is a generator
    loop {
        for p_i in s_factors{
            if mod_pow(g, s/p_i, p) == 1{
                g += 1;
                break;
            }
        }
        return g;
    }
}