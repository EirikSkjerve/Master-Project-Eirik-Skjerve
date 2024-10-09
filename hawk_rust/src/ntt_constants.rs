
use prime_factorization::Factorization;
use crate::utils::{adjoint, mod_pow, modulo};

use std::collections::HashMap;

static mut ZETA_TABLE: Option<HashMap<u128, (Vec<u128>, Vec<u128>)>> = None;

pub fn res_z() {
    unsafe{
        ZETA_TABLE = None;
    }
}

fn z_map() -> &'static mut HashMap<u128, (Vec<u128>, Vec<u128>)> {

    unsafe {

        if ZETA_TABLE.is_none() {
            ZETA_TABLE = Some(HashMap::new());
        }
        return ZETA_TABLE.as_mut().unwrap();
    }

}

fn z_add(key: u128, value: (Vec<u128>, Vec<u128>)) {
    z_map().insert(key,value);
}

fn z_get(key: u128) -> Option<(Vec<u128>, Vec<u128>)> {
    unsafe {
        ZETA_TABLE
            .as_ref()
            .and_then(|map| map.get(&key).cloned())
    }
}

// given prime p and order n, compute 2n-th root of unity mod p
pub fn get_roots(p: u128, n: u128) -> (Vec<u128>, Vec<u128>) {

    if let Some(tbl) = z_get(p) { 
        return tbl;
    } else {
        let mut g0 = primitive_root(p);
        let b = (p - 1) / (2 * n);
        g0 = mod_pow(g0, b, p);

        let zetas = compute_zetas(g0, p, n);
        let mut izetas: Vec<u128> = Vec::new();

        for z in zetas.iter() {
            izetas.push(mod_pow(*z, p - 2, p));
        }
        z_add(p, (zetas, izetas));
    }

    return get_roots(p, n);


}

// comptute zeta values in a specific order for usage in ntt/intt functions
fn compute_zetas(root: u128, p: u128, n: u128) -> Vec<u128> {
    let mut zetas: Vec<u128> = Vec::new();
    let log_n = (n as f64).log2() as u128;

    for u in 0..n {
        let brv = brv(u, log_n);
        zetas.push(mod_pow(root, brv, p));
    }
    return zetas;
}

// return the primitive root of prime p
// using https://www.geeksforgeeks.org/primitive-root-of-a-prime-number-n-modulo-n/
pub fn primitive_root(p: u128) -> u128 {
    // let mut g: u128 = 2;

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
    return 0;
}

// implementation of bit reversal of an integer
pub fn brv(x: u128, log_n: u128) -> u128 {
    let mut r = 0;
    for i in 0..log_n {
        r |= ((x >> i) & 1) << (log_n - 1 - i);
    }
    return r;
}
