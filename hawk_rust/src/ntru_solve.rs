extern crate num_bigint;
extern crate num_traits;

use num_bigint::{BigInt, ToBigInt};
use num_traits::{One, Signed, Zero};

use crate::fft::{add_fft, adj_fft, div_fft, fft, ifft, mul_fft};
use crate::utils::{bigint_to_f64_vec, bigint_to_i64_vec, bigint_vec};
use num_complex::Complex;

// this file contains code for the NTRU-solve algorithm, and its helper-functions
// Because the size of intermediate values in NTRU-solve grows very big, I use the BigInt crate
// for all calculations.
pub fn xgcd(mut a: BigInt, mut b: BigInt) -> (BigInt, BigInt, BigInt) {
    /*
     * Implements extended Euclidean algorithm to find Bezout's coefficients
     * Inputs: integers a and b
     * Outputs: gcd(a,b) and s, t such that a*s + b*t = gcd(a,b)
     */

    // Swap the order if a is less than b
    if a < b {
        let (gcd, t, s) = xgcd(b, a);
        return (gcd, s, t);
    }

    // Initialize coefficients
    let mut s0 = BigInt::one();
    let mut s1 = BigInt::zero();
    let mut t0 = BigInt::zero();
    let mut t1 = BigInt::one();

    // Perform the algorithm
    while !b.is_zero() {
        let q = &a / &b;
        let r = a % &b;

        a = b;
        b = r;

        // Update s and t using in-place operations without clones
        let new_s = &s0 - &q * &s1;
        s0 = s1;
        s1 = new_s;

        let new_t = &t0 - &q * &t1;
        t0 = t1;
        t1 = new_t;
    }

    (a, s0, t0)
}

pub fn karatsuba(a: Vec<BigInt>, b: Vec<BigInt>, n: usize) -> Vec<BigInt> {
    //
    // implementation of the karatsuba algorithm for fast polynomial multiplication
    //

    // base case
    if n == 1 {
        return vec![&a[0] * &b[0], BigInt::ZERO];
    }

    // split polynomials
    let m = n / 2;
    let a0 = a[0..m].to_vec();
    let a1 = a[m..n].to_vec();
    let b0 = b[0..m].to_vec();
    let b1 = b[m..n].to_vec();

    // calculate new factors
    let mut ax = vec![BigInt::ZERO; m];
    let mut bx = vec![BigInt::ZERO; m];
    for i in 0..m {
        ax[i] = &a0[i] + &a1[i];
        bx[i] = &b0[i] + &b1[i];
    }

    // recursive steps
    let c0 = karatsuba(a0, b0, m);
    let c1 = karatsuba(a1, b1, m);
    let mut c2 = karatsuba(ax, bx, m);

    for i in 0..n {
        c2[i] -= &c0[i] + &c1[i];
    }

    // empty vector for keeping the final result
    let mut c = vec![BigInt::ZERO; 2 * n];

    // join terms
    for i in 0..n {
        c[i] += &c0[i];
        c[i + n] += &c1[i];
        c[i + m] += &c2[i];
    }

    c
}

pub fn karamul(a: Vec<BigInt>, b: Vec<BigInt>) -> Vec<BigInt> {
    //
    // perform karatsuba polynomial multiplication with modulo X^n + 1
    //

    // get degree of polynomials
    let n = a.len();

    // perform the multiplication
    let c = karatsuba(a.to_vec(), b.to_vec(), n);

    // reduce result mod X^n + 1
    let mut c_reduced = vec![BigInt::ZERO; n];
    for i in 0..n {
        c_reduced[i] = &c[i] - &c[i + n]
    }
    c_reduced
}

pub fn karamul_sq(a: Vec<BigInt>) -> Vec<BigInt> {
    // squares an element using karatsuba multiplication
    karamul(a.clone(), a.clone())
}

pub fn field_norm(a: Vec<BigInt>) -> Vec<BigInt> {
    //
    // Projects an element from Z[x]/x^n +1 to Z[x]/x^(n/2) + 1
    // by computing n(a) = (af_o)^2 - X(f_e)^2
    // where f_o = f_1 + f_3X + ...
    // and   f_e = f_0 + f_2X + ...
    //

    let m = a.len() / 2;

    // create vectors of size n/2 to keep even and odd coefficients
    let mut a_even: Vec<BigInt> = Vec::with_capacity(m);
    let mut a_odd: Vec<BigInt> = Vec::with_capacity(m);

    // split a into even and odd coefficients
    for i in 0..a.len() {
        if i % 2 == 0 {
            a_even.push(a[i].clone());
        }
        if i % 2 == 1 {
            a_odd.push(a[i].clone());
        }
    }

    // square the polynomials by using karatsuba multiplication
    let a_even_squared = karamul_sq(a_even);
    let a_odd_squared = karamul_sq(a_odd);

    // compute the final results according to the formula
    let mut res = a_even_squared.clone();

    for i in 0..m - 1 {
        res[i + 1] -= &a_odd_squared[i];
    }

    res[0] += &a_odd_squared[m - 1];

    res
}

pub fn lift(a: Vec<BigInt>) -> Vec<BigInt> {
    //
    // compute a lifted version of polynomial a
    // which is equal to a(x^2)
    //

    let n = a.len();

    let mut res: Vec<BigInt> = vec![BigInt::ZERO; 2 * n];
    for i in 0..n {
        res[2 * i] = a[i].clone();
    }
    res
}

pub fn galois_conjugate(a: Vec<BigInt>) -> Vec<BigInt> {
    //
    // calculate galois conjugate of polynomial a
    // which is equal to a(-x)
    //

    let neg_one = -1.to_bigint().unwrap();
    let res: Vec<BigInt> = (0..a.len())
        .map(|i| neg_one.pow(i as u32) * &a[i])
        .collect();
    res
}

pub fn ntrusolve(f: &Vec<i64>, g: &Vec<i64>) -> Option<(Vec<i64>, Vec<i64>)> {
    // solve the NTRU-equation
    // Given f and g, return F and G such that fG - gF = 1,
    // or None if no solution could be found
    //
    // see Prest Pornin 2019 paper
    // https://eprint.iacr.org/2019/015
    //

    // convert polynomials to have coefficients of bigint type
    let f = bigint_vec(f);
    let g = bigint_vec(g);

    // get solution from recursive function
    if let Some((bigf, bigg)) = ntrusolve_inner(&f, &g) {
        // if solution is found, return solution with coefficients as i64
        return Some((bigint_to_i64_vec(bigf), bigint_to_i64_vec(bigg)));
    }
    // return None if no solution were found
    None
}

fn ntrusolve_inner(f: &Vec<BigInt>, g: &Vec<BigInt>) -> Option<(Vec<BigInt>, Vec<BigInt>)> {
    //
    // recursive ntrusolve algorithm
    //

    // get degree
    let n = f.len();

    // innermost recursion is solving the ntru-equation over integers
    // i.e. doing an extended euclidean algorithm
    if n == 1 {
        let (d, u, v) = xgcd(f[0].clone(), g[0].clone());

        // if gcd of f and g at innermost level is not one, we can't find solution to
        // ntru-equation
        if d != BigInt::one() {
            return None;
        } else {
            // otherwise, output from extended euclidean algorithm is solution to ntru-equation
            return Some((vec![-v], vec![u]));
        }
    }

    // project f and g on lower field
    // degree is halved
    let fp = field_norm(f.clone());
    let gp = field_norm(g.clone());

    // solve ntru-equation for halved degree by recursive call
    if let Some((bigfp, biggp)) = ntrusolve_inner(&fp, &gp) {
        // reconstruct solution to ntru-equation for this degree
        // by multiplying F'(x^2) * g(-x)
        // and            G'(x^2) * f(-x)
        let bigf = karamul(lift(bigfp), galois_conjugate(g.clone()));
        let bigg = karamul(lift(biggp), galois_conjugate(f.clone()));

        // reduce the coefficients in the solution polynomials F and G
        let (bigf, bigg) = reduce(f, g, bigf, bigg);

        // return the solution
        return Some((bigf, bigg));
    } else {
        // return None if no solution was found
        return None;
    }
}

pub fn bitsize(a: BigInt) -> u32 {
    //
    // compute (approximate) bitsize of an integer a,
    // up to precision of 8
    //

    // clone a mutable copy of the absolute value of input value
    let mut val = a.abs().clone();
    // create mutable temporary variable
    let mut res: u32 = 0;
    // count bits in intervals of 8
    while val > BigInt::ZERO {
        res += 8;
        val >>= 8;
    }
    res
}

pub fn adjust_fft(
    f: Vec<BigInt>,
    g: Vec<BigInt>,
    size: u32,
) -> (Vec<Complex<f64>>, Vec<Complex<f64>>) {
    //
    // adjusting function for reducing polynomials and converting them to fft representation
    //

    let n = f.len();

    // empty vectors for result
    let mut f_adjust: Vec<BigInt> = Vec::with_capacity(n);
    let mut g_adjust: Vec<BigInt> = Vec::with_capacity(n);

    // each element of new f and g are reduced depending on input parameter "size"
    for i in 0..n {
        f_adjust.push(&f[i] >> (size - 53));
        g_adjust.push(&g[i] >> (size - 53));
    }
    // at this point, f_adjust and g_adjust are small enough to cast into floats
    // compute and return the fft-representation
    let fa_fft = fft(&bigint_to_f64_vec(f_adjust));
    let ga_fft = fft(&bigint_to_f64_vec(g_adjust));

    (fa_fft, ga_fft)
}

pub fn calculate_size(f: Vec<BigInt>, g: Vec<BigInt>) -> u32 {
    // return max of
    //              53
    //              bitsize(max(|f|))
    //              bitsize(max(|g|))
    let f_min = f.iter().min().unwrap().clone();
    let f_max = f.iter().max().unwrap().clone();
    let g_min = g.iter().min().unwrap().clone();
    let g_max = g.iter().max().unwrap().clone();

    let size = [
        53,
        bitsize(f_min),
        bitsize(f_max),
        bitsize(g_min),
        bitsize(g_max),
    ]
    .iter()
    .max()
    .unwrap()
    .clone();

    size
}

pub fn reduce(
    f: &Vec<BigInt>,
    g: &Vec<BigInt>,
    bigf: Vec<BigInt>,
    bigg: Vec<BigInt>,
) -> (Vec<BigInt>, Vec<BigInt>) {
    // since input is vectors of BigInt type, they are not floats. To perform fft-calculations on
    // them, we need to extract the high bits of f, g, F, and G, convert the vectors of the
    // high bits to Vec<f64>, compute a scaling factor k that
    // fits in in i32 type (accounts for sign) using fft.
    // After this, we compute bigf -= kf, bigg -= kg, both as Vec<BigInt>, and return the result.
    // return (f, g, bigf, bigg), with bigf, bigg as reduced;

    let mut bigf_mut = bigf.clone();
    let mut bigg_mut = bigg.clone();

    // find size that determines how much to reduce polynomials with
    let size = calculate_size(f.clone(), g.clone());

    // fft of reduced f and g
    let (fa_fft, ga_fft) = adjust_fft(f.clone(), g.clone(), size);

    loop {
        let size_inner = calculate_size(bigf_mut.clone(), bigg_mut.clone());
        // breaking condition
        if size_inner < size {
            break;
        }

        // fft of reduced f and g
        let (bigfa_fft, bigga_fft) = adjust_fft(bigf_mut.clone(), bigg_mut.clone(), size_inner);

        // calculate ff* + gg*
        let den_fft = add_fft(
            &mul_fft(&fa_fft, &adj_fft(&fa_fft)),
            &mul_fft(&ga_fft, &adj_fft(&ga_fft)),
        );

        // calculate bigff* + biggg*
        let num_fft = add_fft(
            &mul_fft(&bigfa_fft, &adj_fft(&fa_fft)),
            &mul_fft(&bigga_fft, &adj_fft(&ga_fft)),
        );

        // calculate factor k, first through fft
        let k_f = ifft(&div_fft(&num_fft, &den_fft));

        // convert to integers with rounding
        let k: Vec<BigInt> = bigint_vec(&(0..k_f.len()).map(|x| k_f[x].round() as i64).collect());

        // check if k is zero
        let mut zero_flag = true;
        for i in 0..k.len() {
            if k[i] != BigInt::ZERO {
                zero_flag = false;
            }
        }
        // if so, break
        if zero_flag {
            break;
        }

        // calculate factors that should be reduced from bigf and bigg
        let fk = karamul(f.clone(), k.clone());
        let gk = karamul(g.clone(), k.clone());

        // reduce bigf and bigg
        for i in 0..f.len() {
            bigf_mut[i] -= fk[i].clone() << (size_inner - size);
            bigg_mut[i] -= gk[i].clone() << (size_inner - size);
        }
    }
    (bigf_mut, bigg_mut)
}
