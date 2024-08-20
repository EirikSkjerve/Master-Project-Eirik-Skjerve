/*

pub fn ntrusolve(f: &Vec<i64>, g: &Vec<i64>) -> (Vec<i64>, Vec<i64>, Vec<i64>, Vec<i64>){


    // return (f, g, f, g);
}

pub fn reduce(f: &Vec<i64>, g: &Vec<i64>, G: &Vec<i64>, G: &Vec<i64>) -> (Vec<i64>, Vec<i64>, Vec<i64>, Vec<i64>){

    // return (f, g, F, G);
}
*/
use crate::utils::modulo;

pub fn xgcd(a_inp: i64, b_inp: i64) -> (i64, i64, i64){
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
    let mut cof: [i64; 4] = [1, 0, 0, 1];
    let mut a = a_inp;
    let mut b = b_inp;
    let mut q = 0;
    let mut temp = 0;

    // perform the algorithm
    while b != 0 {

        // rounded division
        q = a/b;

        // calculates the gcd
        temp = b;
        b = modulo(a, b);
        a = temp;

        // calculates the coefficients
        temp = cof[1];
        cof[1] = cof[0] - q*cof[1];
        cof[0] = temp;

        temp = cof[3];
        cof[3] = cof[2] - q*cof[3];
        cof[2] = temp;
    }

    return (a, cof[0], cof[2]);
}

// implementation of the karatsuba algorithm for fast integer/polynomial multiplication
pub fn karatsuba(a: &Vec<i64>, b: &Vec<i64>, n: usize) -> Vec<i64> {

    // base case
    if n==1{
        return vec![a[0]*b[0], 0];
    }

    let m = n/2;
    let a0 = a[0..m].to_vec();
    let a1 = a[m..n as usize].to_vec();
    let b0 = b[0..m].to_vec();
    let b1 = b[m..n as usize].to_vec();

    let mut ax = vec![0;m];
    let mut bx = vec![0;m];
    for i in 0..m {
        ax[i] = a0[i] + a1[i];
        bx[i] = b0[i] + b1[i];
    }

    let c0 = karatsuba(&a0, &b0, m);
    let c1 = karatsuba(&a1, &b1, m);

    let mut c2 = karatsuba(&ax, &bx, m);


    for i in 0..n{
        c2[i] -= c0[i] + c1[i];
    }

    let mut c = vec![0; 2*n];

    for i in 0..n {
        c[i] += c0[i];
        c[i+n] += c1[i];
        c[i+m] += c2[i];
    }

    return c;
}

pub fn karamul(a: &Vec<i64>, b: &Vec<i64>) -> Vec<i64> {
    let n = a.len();
    let c = karatsuba(a, b, n);
    let mut c_reduced = vec![0; n];

    for i in 0..n {
        c_reduced[i] = c[i] - c[i+n]
    }
    return c_reduced;
}

pub fn field_norm(a: &Vec<i64>) -> Vec<i64> {
    /*
     * Projects an element from Q[x]/x^n +1 to Q[x]/x^(n/2) + 1
     */
    let m = a.len()/2;

    // split a into even and odd coefficients
    let mut a_even: Vec<i64> = Vec::new();
    let mut a_odd: Vec<i64> = Vec::new();

    for i in 0..a.len() {
        if i%2 == 0{
            a_even.push(a[i]);
        }
        if i%2 == 1 {
            a_odd.push(a[i]);
        }
    }

    // square polynomials
    let a_even_squared = karamul(&a_even, &a_even);
    let a_odd_squared = karamul(&a_odd, &a_odd);

    let mut res = a_even_squared.clone();
    println!("{:?}", a_even_squared);


    for i in 0..m-1 {
        //res[i+1] = a_even_squared[i+1] - a_odd_squared[i];
        res[i+1] -= a_odd_squared[i];
    }

    res[0] += a_odd_squared[m-1];
    return res;
}

