// This is the Number Theoretic Transform module
// see https://eprint.iacr.org/2024/585.pdf
// and section 4.1.1 of HAWK spec paper

use crate::utils::{adjoint, mod_pow, modulo};
use crate::ntt_constants::get_roots;
use std::time::Instant;

// ntt of a polynomial
pub fn ntt(f: Vec<i64>, p: u32) -> Vec<i64> {
    let mut ntt_f = f.clone();

    let n = f.len();
    let mut l = n / 2;
    let mut k = 1;


    let m = Instant::now();
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

pub fn nttadj(f: &Vec<i64>, p: u32) -> Vec<i64> {
    let ui = intt(f.to_vec(), p);
    return ntt(adjoint(&ui), p);
}



