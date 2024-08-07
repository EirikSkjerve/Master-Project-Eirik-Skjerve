// This is the Number Theoretic Transform module
// see https://eprint.iacr.org/2024/585.pdf
// and section 4.1.1 of HAWK spec paper

// use num::One;
use crate::utils::{mod_pow, modulo};
use num::traits::{FromPrimitive, Num, PrimInt};
use prime_factorization::Factorization;

// ntt of a polynomial
pub fn ntt(f: Vec<i64>, p: u32) -> Vec<i64> {
    let mut ntt_f = f.clone();

    let n = f.len();
    let mut l = n / 2;
    let mut k = 1;

    let zetas = get_roots(p as u128, n as u128).0;

    let q = p as i64;

    while l > 0 {
        for s in (0..n).step_by(2 * l) {
            let zeta = zetas[k] as i64;
            k += 1;
            for j in s..s + l {
                let t = modulo(ntt_f[j + l] * zeta, q);
                ntt_f[j + l] = modulo(ntt_f[j] - t, q);
                ntt_f[j] = modulo(ntt_f[j] + t, q);
            }
        }
        l /= 2;
    }

    return ntt_f;
}

// inverse ntt of a polynomial
pub fn intt(f: Vec<i64>, p: u32) -> Vec<i64> {
    let mut intt_f = f.clone();

    let n = f.len();
    let mut l = 1;
    let mut k = n - 1;

    let izetas = get_roots(p as u128, n as u128).1;

    let q = p as i64;

    while l < n {
        for s in (0..n).step_by(2 * l).rev() {
            let izeta = izetas[k] as i64;
            k -= 1;
            for j in s..s + l {
                let t = intt_f[j];
                intt_f[j] = modulo(t + intt_f[j + l], q);
                intt_f[j + l] = modulo(t - intt_f[j + l], q);
                intt_f[j + l] = modulo(intt_f[j + l] * izeta, q);
            }
        }
        l *= 2;
    }

    let ideg = mod_pow(n as i64, q - 2, q);
    for i in 0..intt_f.len() {
        intt_f[i] = modulo(intt_f[i] * ideg, q);
    }

    return intt_f;
}

// implementation of bit reversal of an integer
pub fn brv(x: u128, log_n: u128) -> u128 {
    // following is HAWK-teams implementation. Might go about it a different way
    // no idea how it works though

    let mut r = 0;
    for i in 0..log_n {
        r |= ((x >> i) & 1) << (log_n - 1 - i);
    }
    return r;
}

// given prime p and order n, compute 2n-th root of unity mod p
pub fn get_roots(p: u128, n: u128) -> (Vec<u128>, Vec<u128>) {
    let mut g0 = primitive_root(p);
    let b = (p - 1) / (2 * n);
    g0 = mod_pow(g0, b, p);

    let zetas = compute_zetas(g0, p, n);
    let mut izetas: Vec<u128> = Vec::new();

    for z in zetas.iter() {
        izetas.push(mod_pow(*z, p - 2, p));
    }

    return (zetas, izetas);
}

// comptute zeta values in a specific order for usage in ntt/intt functions
fn compute_zetas(root: u128, p: u128, n: u128) -> Vec<u128> {
    let mut zetas: Vec<u128> = Vec::new();
    let log_n = (n as f64).log2() as u128;

    for u in 0..n {
        let brv = brv(u, log_n);
        //println!("brv: {}", brv);
        let res = mod_pow(root, brv, p);
        //println!("res: {}", res);
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
    let mut flag = false;
    for g in 2..p {
        flag = false;
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
