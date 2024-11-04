// This is the Number Theoretic Transform module
// see https://eprint.iacr.org/2024/585.pdf
// and section 4.1.1 of HAWK spec paper

use crate::ntt_constants::get_roots;
use crate::utils::{adjoint, mod_pow, modulo};

pub fn ntt(f: Vec<i64>, p: i64) -> Vec<i64> {
    //
    // compute the Number Theoretic Transform (NTT) representation of a polynomial f
    // requires some prime p
    //

    // make a copy of f to store the new coefficients in
    let mut ntt_f = f.clone();

    let n = f.len();

    // counter variables for the algorithm
    let mut l = n / 2;
    let mut k = 1;

    // compute zetas and izetas for input p and n
    // these will only be computed the first time get_roots() is called;
    // otherwise they will be ready in pre-computed hash-table
    let zetas = get_roots(p as u128, n as u128).0;

    // perform the algorithm as described in section 4.1.1 of Hawk spec paper
    while l > 0 {
        for s in (0..n).step_by(2 * l) {
            let zeta = zetas[k] as i64;
            k += 1;
            for j in s..s + l {
                let t = modulo(ntt_f[j + l] * zeta, p);
                ntt_f[j + l] = modulo(ntt_f[j] - t, p);
                ntt_f[j] = modulo(ntt_f[j] + t, p);
            }
        }
        l /= 2;
    }

    ntt_f
}

pub fn intt(f: Vec<i64>, p: i64) -> Vec<i64> {
    //
    // compute the inverse Number Theoretic Transform (iNTT) of a polynomial in NTT representation
    // requires some prime p
    //
    // This is essentially the inverse function of ntt(f, p), i.e. intt(ntt(f, p), p) = f
    //

    // make a copy of f to store new coefficients in
    let mut intt_f = f.clone();

    let n = f.len();

    // counter variables for the algorithm
    let mut l = 1;
    let mut k = n - 1;

    // see comment in ntt()
    let izetas = get_roots(p as u128, n as u128).1;

    // perform the algorithm as described in section 4.1.1 of Hawk spec paper
    while l < n {
        for s in (0..n).step_by(2 * l).rev() {
            let izeta = izetas[k] as i64;
            k -= 1;
            for j in s..s + l {
                let t = intt_f[j];
                intt_f[j] = modulo(t + intt_f[j + l], p);
                intt_f[j + l] = modulo(t - intt_f[j + l], p);
                intt_f[j + l] = modulo(intt_f[j + l] * izeta, p);
            }
        }
        l *= 2;
    }

    // final step to make coefficients correct
    let ideg = mod_pow(n as i64, p - 2, p);
    for i in 0..intt_f.len() {
        intt_f[i] = modulo(intt_f[i] * ideg, p);
    }

    intt_f
}

pub fn nttadj(f: &Vec<i64>, p: i64) -> Vec<i64> {
    //
    // compute the NTT representation of hermitian adjoint of f
    // used for computing public key Q in keygen
    //

    let ui = intt(f.to_vec(), p);
    ntt(adjoint(&ui), p)
}
