extern crate num_bigint;
extern crate num_traits;

use num_bigint::{BigInt, BigUint, ToBigInt, ToBigUint};
use num_traits::{One, Signed, ToPrimitive, Zero};

use crate::fft::{add_fft, adj_fft, div_fft, fft, ifft, mul_fft, sub_fft};
use crate::utils::bigint_vec;
use num_complex::Complex;

// this file contains code for the NTRU-solve algorithm, and its helper-functions
// Because the size of intermediate values in NTRU-solve grows very big, I use the BigInt crate
// for all calculations.

pub fn xgcd(a_inp: BigInt, b_inp: BigInt) -> (BigInt, BigInt, BigInt) {
    /*
     * Implements extended euclidean algorithm to find Bezout's coefficients
     * Inputs: integers a and b
     * Outputs: gcd(a,b) and s, t such that a*s + b*t = gcd(a,b)
     */

    // swap the order if a is less than b
    if a_inp < b_inp {
        let res = xgcd(b_inp, a_inp);
        return (res.0, res.2, res.1);
    }

    // initialize variables
    let mut cof: [BigInt; 4] = [BigInt::one(), BigInt::ZERO, BigInt::ZERO, BigInt::one()];
    let mut a = a_inp;
    let mut b = b_inp;
    let mut q = BigInt::ZERO;
    let mut temp = BigInt::ZERO;

    // perform the algorithm
    while b != BigInt::ZERO {
        // rounded division
        q = &a / &b;

        // calculates the gcd
        temp = b.clone();
        b = a % &b;
        a = temp;

        // calculates the coefficients
        temp = cof[1].clone();
        cof[1] = cof[0].clone() - &q * cof[1].clone();
        cof[0] = temp;

        temp = cof[3].clone();
        cof[3] = cof[2].clone() - &q * cof[3].clone();
        cof[2] = temp;
    }

    return (a, cof[0].clone(), cof[2].clone());
}

// implementation of the karatsuba algorithm for fast integer/polynomial multiplication
pub fn karatsuba(a: Vec<BigInt>, b: Vec<BigInt>, n: usize) -> Vec<BigInt> {
    // base case
    if n == 1 {
        return vec![&a[0] * &b[0], BigInt::ZERO];
    }

    let m = n / 2;
    let a0 = a[0..m].to_vec();
    let a1 = a[m..n as usize].to_vec();
    let b0 = b[0..m].to_vec();
    let b1 = b[m..n as usize].to_vec();

    let mut ax = vec![BigInt::ZERO; m];
    let mut bx = vec![BigInt::ZERO; m];
    for i in 0..m {
        ax[i] = &a0[i] + &a1[i];
        bx[i] = &b0[i] + &b1[i];
    }

    let c0 = karatsuba(a0, b0, m);
    let c1 = karatsuba(a1, b1, m);

    let mut c2 = karatsuba(ax, bx, m);

    for i in 0..n {
        c2[i] -= &c0[i] + &c1[i];
    }

    let mut c = vec![BigInt::ZERO; 2 * n];

    for i in 0..n {
        c[i] += &c0[i];
        c[i + n] += &c1[i];
        c[i + m] += &c2[i];
    }

    return c;
}

pub fn karamul(a: Vec<BigInt>, b: Vec<BigInt>) -> Vec<BigInt> {
    let n = a.len();
    let c = karatsuba(a.to_vec(), b.to_vec(), n);
    let mut c_reduced = vec![BigInt::ZERO; n];

    for i in 0..n {
        c_reduced[i] = &c[i] - &c[i + n]
    }
    return c_reduced;
}

// squares an element using karatsuba multiplication
pub fn karamul_sq(a: Vec<BigInt>) -> Vec<BigInt> {
    return karamul(a.clone(), a.clone());
}

pub fn field_norm(a: Vec<BigInt>) -> Vec<BigInt> {
    /*
     * Projects an element from Q[x]/x^n +1 to Q[x]/x^(n/2) + 1
     */
    let m = a.len() / 2;

    // split a into even and odd coefficients
    let mut a_even: Vec<BigInt> = Vec::new();
    let mut a_odd: Vec<BigInt> = Vec::new();

    for i in 0..a.len() {
        if i % 2 == 0 {
            a_even.push(a[i].clone());
        }
        if i % 2 == 1 {
            a_odd.push(a[i].clone());
        }
    }

    // square polynomials
    let a_even_squared = karamul_sq(a_even);
    let a_odd_squared = karamul_sq(a_odd);

    let mut res = a_even_squared.clone();

    for i in 0..m - 1 {
        res[i + 1] -= &a_odd_squared[i];
    }

    res[0] += &a_odd_squared[m - 1];
    return res;
}

pub fn lift(a: Vec<BigInt>) -> Vec<BigInt> {
    let n = a.len();

    let mut res: Vec<BigInt> = vec![BigInt::ZERO; 2 * n];
    for i in 0..n {
        res[2 * i] = a[i].clone();
    }
    return res;
}

pub fn galois_conjugate(a: Vec<BigInt>) -> Vec<BigInt> {
    let n = a.len();
    let neg_one = -1.to_bigint().unwrap();
    let res: Vec<BigInt> = (0..n).map(|i| neg_one.pow(i as u32) * &a[i]).collect();
    return res;
}

pub fn ntrusolve(f: Vec<BigInt>, g: Vec<BigInt>) -> (Vec<BigInt>, Vec<BigInt>) {
    let n = f.len();

    if n == 1 {
        let (d, u, v) = xgcd(f[0].clone(), g[0].clone());

        if d != BigInt::one() {
            println!("gcd({}, {}) = {}, aborting", f[0], g[0], d);
            // this should throw an error or return false
            return (vec![BigInt::ZERO], vec![BigInt::ZERO]);
        } else {
            // in theory, this return -q*v, q*u, but in HAWK, q=1
            return (vec![-v], vec![u]);
        }
    }

    let fp = field_norm(f.clone());
    let gp = field_norm(g.clone());

    let (Fp, Gp) = ntrusolve(fp, gp);

    let mut F = karamul(lift(Fp), galois_conjugate(g.clone()));
    let mut G = karamul(lift(Gp), galois_conjugate(f.clone()));

    // println!("before: F: {:?}, G:{:?} \n", F, G);
    let (F, G) = reduce(f, g, F, G);
    // println!("\nafter: F: {:?}, G:{:?}", F, G);
    return (F, G);
}

pub fn bitsize(a: BigInt) -> u32 {
    // clone a mutable copy of the absolute value of input value
    let mut val = a.abs().clone();

    // create mutable temporary variable
    let mut res: u32 = 0;
    // count bits in intervals of 8
    while val > BigInt::ZERO {
        res += 8;
        val >>= 8;
    }
    // try and convert the temporary value into an u32
    if let Some(res_u32) = res.to_u32() {
        return res_u32;
    } else {
        println!("res_u32 variable is too big");
        return 0;
    }
}

pub fn bigint_to_f64_vec(a: Vec<BigInt>) -> Vec<f64> {
    let n = a.len();
    let mut res: Vec<f64> = Vec::with_capacity(n);

    for i in 0..n {
        if let Some(res_i) = &a[i].to_f64() {
            res.push(*res_i);
        } else {
            println!("Could not convert to float");
        }
    }

    return res;
}

pub fn adjust_fft(
    f: Vec<BigInt>,
    g: Vec<BigInt>,
    size: u32,
) -> (Vec<Complex<f64>>, Vec<Complex<f64>>) {
    let n = f.len();
    // set the minimum threshold bitsize to
    let mut f_adjust: Vec<BigInt> = Vec::with_capacity(n);
    let mut g_adjust: Vec<BigInt> = Vec::with_capacity(n);

    for i in 0..n {
        f_adjust.push(&f[i] >> (size - 53));
        g_adjust.push(&g[i] >> (size - 53));
    }
    // at this point, f_adjust and g_adjust should be small enough to cast into floats
    let fa_fft = fft(&bigint_to_f64_vec(f_adjust));
    let ga_fft = fft(&bigint_to_f64_vec(g_adjust));

    return (fa_fft, ga_fft);
}

pub fn calculate_size(f: Vec<BigInt>, g: Vec<BigInt>) -> u32 {
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

    return size;
}

pub fn reduce(
    f: Vec<BigInt>,
    g: Vec<BigInt>,
    F: Vec<BigInt>,
    G: Vec<BigInt>,
) -> (Vec<BigInt>, Vec<BigInt>) {
    // Pornin, Prest 2019's method, also used in HAWK's implementation, is the following:
    // since input is vectors of BigInt type, they are not floats. To perform fft-calculations on
    // them, we need to extract the high bits of f, g, F, and G, convert the vectors of the
    // high bits to Vec<f32> or Vec<f64>, compute a scaling factor k that
    // fits in in i32 type (accounts for sign) using fft.
    // After this, we compute F -= kf, G -= kg, both as Vec<BigInt>, and return the result.
    // return (f, g, F, G);

    let mut F_mut = F.clone();
    let mut G_mut = G.clone();

    let size = calculate_size(f.clone(), g.clone());
    let (fa_fft, ga_fft) = adjust_fft(f.clone(), g.clone(), size);

    // correct up to here

    loop {
        let size_inner = calculate_size(F_mut.clone(), G_mut.clone());
        if size_inner < size {
            break;
        }

        let (Fa_fft, Ga_fft) = adjust_fft(F_mut.clone(), G_mut.clone(), size_inner);

        // calculate ff* + gg*
        let den_fft = add_fft(
            &mul_fft(&fa_fft, &adj_fft(&fa_fft)),
            &mul_fft(&ga_fft, &adj_fft(&ga_fft)),
        );
        // calculate Ff* + Gg*
        let num_fft = add_fft(
            &mul_fft(&Fa_fft, &adj_fft(&fa_fft)),
            &mul_fft(&Ga_fft, &adj_fft(&ga_fft)),
        );

        // calculate factor k, first through fft
        let k_f = ifft(&div_fft(&num_fft, &den_fft));

        // convert to integers
        let k: Vec<BigInt> = bigint_vec((0..k_f.len()).map(|x| k_f[x].round() as i64).collect());

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

        // calculate factors that should be reduced from F and G
        let fk = karamul(f.clone(), k.clone());
        let gk = karamul(g.clone(), k.clone());

        // reduce F and G
        for i in 0..f.len() {
            F_mut[i] -= fk[i].clone() << (size_inner - size);
            G_mut[i] -= gk[i].clone() << (size_inner - size);
        }
    }
    return (F_mut, G_mut);
}
