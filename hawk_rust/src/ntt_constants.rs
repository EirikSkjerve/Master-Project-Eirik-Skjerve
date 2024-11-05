use crate::utils::mod_pow;
use prime_factorization::Factorization;

use std::collections::HashMap;

// this file contains functionality for computing zetas and izetas used in the ntt module

// initialize empty hashmap
static mut ZETA_TABLE: Option<HashMap<u128, (Vec<u128>, Vec<u128>)>> = None;

// unsure if this is needed
// pub fn res_z() {
//     unsafe {
//         ZETA_TABLE = None;
//     }
// }

fn z_map() -> &'static mut HashMap<u128, (Vec<u128>, Vec<u128>)> {
    //
    // function for getting a reference to the hashmap keeping the roots
    //
    unsafe {
        // if hashmap is not yet initialized, initialize it
        if ZETA_TABLE.is_none() {
            ZETA_TABLE = Some(HashMap::new());
        }

        // return a reference to the hashmap
        return ZETA_TABLE.as_mut().unwrap();
    }
}

fn z_add(key: u128, value: (Vec<u128>, Vec<u128>)) {
    // insert an entry into the hashmap
    z_map().insert(key, value);
}

fn z_get(key: u128) -> Option<(Vec<u128>, Vec<u128>)> {
    // get a value from the hashmap
    unsafe { ZETA_TABLE.as_ref().and_then(|map| map.get(&key).cloned()) }
}

pub fn get_roots(p: u128, n: u128) -> (Vec<u128>, Vec<u128>) {
    //
    // given prime p and order n, compute 2n-th roots of unity mod p
    // the first time this function is called with p, the roots will be computed from scratch
    // otherwise this will only return a precomputed table of roots
    //
    if let Some(tbl) = z_get(p) {
        return tbl;
    } else {
        // compute primitive root of prime p
        // and root of unity
        let mut g0 = primitive_root(p);
        let b = (p - 1) / (2 * n);
        g0 = mod_pow(g0, b, p);

        // compute zetas
        let zetas = compute_zetas(g0, p, n);

        // compute izetas given zetas
        let mut izetas: Vec<u128> = Vec::new();

        for z in zetas.iter() {
            izetas.push(mod_pow(*z, p - 2, p));
        }

        // add entry to hashmap
        z_add(p, (zetas, izetas));
    }

    // return a recursive call after creating the hashmap
    get_roots(p, n)
}

fn compute_zetas(root: u128, p: u128, n: u128) -> Vec<u128> {
    // comptute zeta values in a specific order for usage in ntt/intt functions
    let mut zetas: Vec<u128> = Vec::new();
    let log_n = (n as f64).log2() as u128;

    for u in 0..n {
        let brv = brv(u, log_n);
        zetas.push(mod_pow(root, brv, p));
    }
    zetas
}

pub fn primitive_root(p: u128) -> u128 {
    // return the primitive root of prime p
    // using https://www.geeksforgeeks.org/primitive-root-of-a-prime-number-n-modulo-n/

    // phi(p) = p-1
    let s = p - 1;

    // compute prime factors of p-1
    // thanks to crate
    let mut s_factors = Factorization::run(s).factors;

    // remove duplicates
    s_factors.dedup();

    // check if g is a generator
    for g in 2..p {
        let mut flag = false;
        for p_i in s_factors.iter() {
            if mod_pow(g, s / p_i, p) == 1 {
                flag = true;
            }
        }
        // return g if flag is not raised <=> g is primitive root
        if !flag {
            return g;
        }
    }
    // return default value so compiler is happy
    0
}

pub fn brv(x: u128, log_n: u128) -> u128 {
    // implementation of bit reversal of an integer
    let mut r = 0;
    for i in 0..log_n {
        r |= ((x >> i) & 1) << (log_n - 1 - i);
    }
    r
}
