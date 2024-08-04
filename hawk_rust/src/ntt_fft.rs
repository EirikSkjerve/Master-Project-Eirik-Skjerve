// This is the Number Theoretic Transform module
// see https://eprint.iacr.org/2024/585.pdf

use prime_factorization::Factorization;
use crate::utils::mod_pow;

// ntt of a polynomial
pub fn ntt(f: Vec<i32>, p: u32) -> Vec<i32> {

    let mut ntt_f = vec![0;f.len()];

    // using *this* algorithm

    return ntt_f;
}

// implementation of bit reversal of an integer
fn brv(n: i32) -> i32 {
    // assert n power of 2
    // assert b = log2(n)
    return 0;
}

// given prime p and order n, compute 2n-th root of unity mod p
fn get_roots(p: i32, n: i32) -> (Vec<i32>, Vec<i32>) {

    let mut g0 = primitive_root(p);
    let b = (p-1)/(2*n);
    g0 = mod_pow(g0, b, p);
    println!("g0: {}", g0);

    let zetas = compute_zetas(g0, p, n);    
    let mut izetas: Vec<i32> = Vec::new();

    for z in zetas.iter() {
        izetas.push(mod_pow(*z, p-2, p));
    }

    return (zetas, izetas);
}

// comptute zeta values in a specific order for usage in ntt/intt functions
fn compute_zetas(root: i32, p: i32, n: i32) -> Vec<i32>{
    let zetas = vec![0; n as usize];

    let log_n = (n as f64).log2() as i32;
    return zetas;
}

// return the primitive root of prime p
// using https://www.geeksforgeeks.org/primitive-root-of-a-prime-number-n-modulo-n/
pub fn primitive_root(p: i32) -> i32 {
    let mut g = 2;

    // phi(p) = p-1
    let s = p-1;

    // compute prime factors of p-1
    let mut s_factors = Factorization::run(s as u64).factors;

    // remove duplicates
    s_factors.dedup();

    // check if g is a generator
    loop {
        for p_i in s_factors{
            if mod_pow(g, (s as usize/p_i as usize) as u64, p as u64) == 1{
                g += 1;
                break;
            }
        }
        println!("Generator for {}: {}", p, g);
        return g as i32;
    }
}
