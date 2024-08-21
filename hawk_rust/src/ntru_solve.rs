extern crate num_bigint;
extern crate num_traits;

use num_bigint::{BigInt, BigUint, ToBigInt, ToBigUint};
use num_traits::{One, Zero, ToPrimitive, Signed};
use crate::utils::modulo;


// this file contains code for the NTRU-solve algorithm, and its helper-functions
// Because the size of intermediate values in NTRU-solve grows very big, I use the BigInt crate
// for all calculations.

pub fn xgcd(a_inp: BigInt, b_inp: BigInt) -> (BigInt, BigInt, BigInt){
    /*
     * Implements extended euclidean algorithm to find Bezout's coefficients
     * Inputs: integers a and b
     * Outputs: gcd(a,b) and s, t such that a*s + b*t = gcd(a,b)
    */

    // swap the order if a is less than b
    if a_inp < b_inp{
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
        q = &a/&b;

        // calculates the gcd
        temp = b.clone();
        b = a % &b;
        a = temp;

        // calculates the coefficients
        temp = cof[1].clone();
        cof[1] = cof[0].clone() - &q*cof[1].clone();
        cof[0] = temp;

        temp = cof[3].clone();
        cof[3] = cof[2].clone() - &q*cof[3].clone();
        cof[2] = temp;
    }

    return (a, cof[0].clone(), cof[2].clone());
}

// implementation of the karatsuba algorithm for fast integer/polynomial multiplication
pub fn karatsuba(a: Vec<BigInt>, b: Vec<BigInt>, n: usize) -> Vec<BigInt> {

    // base case
    if n==1{
        return vec![&a[0]*&b[0], BigInt::ZERO];
    }

    let m = n/2;
    let a0 = a[0..m].to_vec();
    let a1 = a[m..n as usize].to_vec();
    let b0 = b[0..m].to_vec();
    let b1 = b[m..n as usize].to_vec();

    let mut ax = vec![BigInt::ZERO;m];
    let mut bx = vec![BigInt::ZERO;m];
    for i in 0..m {
        ax[i] = &a0[i] + &a1[i];
        bx[i] = &b0[i] + &b1[i];
    }

    let c0 = karatsuba(a0, b0, m);
    let c1 = karatsuba(a1, b1, m);

    let mut c2 = karatsuba(ax, bx, m);


    for i in 0..n{
        c2[i] -= &c0[i] + &c1[i];
    }

    let mut c = vec![BigInt::ZERO; 2*n];

    for i in 0..n {
        c[i] += &c0[i];
        c[i+n] += &c1[i];
        c[i+m] += &c2[i];
    }

    return c;
}

pub fn karamul(a: Vec<BigInt>, b: Vec<BigInt>) -> Vec<BigInt> {
    let n = a.len();
    let c = karatsuba(a.to_vec(), b.to_vec(), n);
    let mut c_reduced = vec![BigInt::ZERO; n];

    for i in 0..n {
        c_reduced[i] = &c[i] - &c[i+n]
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
    let m = a.len()/2;

    // split a into even and odd coefficients
    let mut a_even: Vec<BigInt> = Vec::new();
    let mut a_odd: Vec<BigInt> = Vec::new();

    for i in 0..a.len() {
        if i%2 == 0{
            a_even.push(a[i].clone());
        }
        if i%2 == 1 {
            a_odd.push(a[i].clone());
        }
    }

    // square polynomials
    let a_even_squared = karamul_sq(a_even);
    let a_odd_squared = karamul_sq(a_odd);

    let mut res = a_even_squared.clone();

    for i in 0..m-1 {
        res[i+1] -= &a_odd_squared[i];
    }

    res[0] += &a_odd_squared[m-1];
    return res;
}

pub fn lift(a: Vec<BigInt>) -> Vec<BigInt> {
    let n = a.len();

    let mut res: Vec<BigInt> = vec![BigInt::ZERO; 2*n];
    for i in 0..n {
        res[2*i] = a[i].clone();
    }
    return res;
}

pub fn galois_conjugate(a: Vec<BigInt>) -> Vec<BigInt> {
    let n = a.len();
    // let res: Vec<BigInt> = vec![BigInt::ZERO; 2*n];
    let neg_one = -1.to_bigint().unwrap();
    let res: Vec<BigInt> = (0..n).map(|i| neg_one.pow(i as u32) * &a[i]).collect();
    return res;
}

pub fn ntrusolve(f: Vec<BigInt>, g: Vec<BigInt>) -> (Vec<BigInt>, Vec<BigInt>){

    let n = f.len();

    if n == 1 {
        let (d, u, v) = xgcd(f[0].clone(), g[0].clone());

        if d != BigInt::one(){
            println!("gcd({}, {}) = {}, aborting",f[0], g[0], d);
            // this should throw an error or return false
            return (vec![BigInt::ZERO], vec![BigInt::ZERO]);
        }
        else {
            // in theory, this return -q*v, q*u, but in HAWK, q=1
            return (vec![-v], vec![u]);
        }
    }

    let fp = field_norm(f.clone());
    let gp = field_norm(g.clone());

    let (Fp, Gp) = ntrusolve(fp, gp);

    let mut F = karamul(lift(Fp), galois_conjugate(g));
    let mut G = karamul(lift(Gp), galois_conjugate(f));
    
    // let (F, G) = reduce(f, g, F, G);

    return (F, G);
}

pub fn bitsize(a: BigInt) -> u32 {
    // clone a mutable copy of the absolute value of input value
    let mut val = a.abs().clone();
   
    // create mutable temporary variable
    let mut res: u32 = 0;
    // count bits in intervals of 8
    while val > BigInt::ZERO{
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

pub fn reduce(f: Vec<BigInt>, g: Vec<BigInt>, F: Vec<BigInt>, G: Vec<BigInt>) -> (Vec<BigInt>, Vec<BigInt>){

    // Pornin, Prest 2019's method, also used in HAWK's implementation, is the following:
    // since input is vectors of BigInt type, they are not floats. To perform fft-calculations on
    // them, we need to extract the high bits of f, g, F, and G, convert the vectors of the 
    // high bits to Vec<f32> or Vec<f64>, compute a scaling factor k that
    // fits in in i32 type (accounts for sign) using fft.
    // After this, we compute F -= kf, G -= kg, both as Vec<BigInt>, and return the result.
    // return (f, g, F, G);
    
    let n = f.len();
    // set the minimum threshold bitsize to 
    let f_min = f.iter().min().unwrap().clone();
    let f_max = f.iter().max().unwrap().clone();
    let g_min = g.iter().min().unwrap().clone();
    let g_max = g.iter().max().unwrap().clone();

    let size = [53, bitsize(f_min), bitsize(f_max), bitsize(g_min), bitsize(g_max)].iter().max().unwrap();
    
    return (F, G);
}
