// This is the Number Theoretic Transform module
// see https://eprint.iacr.org/2024/585.pdf

use crate::utils::fbe;
use prime_factorization::Factorization;

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
fn get_roots(p: i32, n: i32) -> Vec<i32>{
    let roots = vec![0; (n as usize)];
    return roots;
}

// comptute zeta values in a specific order for usage in ntt/intt functions
fn compute_zetas(root: i32, p: i32, n: i32) -> Vec<i32>{
    let zetas = vec![0; (n as usize)];
    return zetas;
}

// return the primitive root of prime p
// using https://www.geeksforgeeks.org/primitive-root-of-a-prime-number-n-modulo-n/
pub fn primitive_root(p: i32) -> i32 {
    let g = 2;
    // phi(p) = p-1
    let s = p-1;
    // compute prime factors of p-1
    let mut p_factors = Factorization::run(s as u64).factors;
    //let p_factors = p_factors.factors.into_iter().unique().collect();
    p_factors.dedup();
    println!("{:?}", p_factors);
    return g;
}
