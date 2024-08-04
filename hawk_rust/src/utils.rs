use num::traits::{Num, FromPrimitive, PrimInt};
use crate::ntt_fft;

pub fn bin(a: u128, x: usize) -> Vec<u8> {
    /*
    Converts an integer to binary representation in a vector of arbitrary size
    Call bin(a, x=0) for default binary size
    */
    let mut res: Vec<u8> = vec![];
    let mut b = a;

    if x == 0 {
        // returns a as binary with *default* size
        return bin(a, (a as f64).log2().ceil() as usize + 1)
    }
    for i in 0..x {
        if b % 2 == 1 {
            res.push(1);
        }
        if b % 2 == 0 {
            res.push(0);
        }
        b /= 2;
    }

    res.reverse();
    return res;
}

pub fn int<T: AsRef<[u8]>>(input: T) -> u128 {
    let a = input.as_ref();
    /*
    Converts an array representing a binary string to its decimal representation
    and returns an 128 bit integer
    */
    let mut res: u128 = 0;
    for i in 0..a.len() {
        if a[i] == 1 {
            res += (2_u128).pow((a.len() - (i + 1)) as u32);
        }
    }
    return res;
}

// implements fast binary exponentiation for computing base^exp mod modulus
// inputs base, exponent and modulus as generic, and returns a u128
pub fn mod_pow<T: PrimInt>(base: T, exp: T, modulus: T) -> T
where 
    T: Num + FromPrimitive,
{
    // convert the inputs to u64 
    let mut base_u64 = base.to_u64().unwrap();
    let mut exp_u64 = exp.to_u64().unwrap();
    let mod_u64 = modulus.to_u64().unwrap();

    // perform the algorithm
    let mut result = 1;
    base_u64 %= mod_u64;
    while exp_u64 > 0 {
        if exp_u64 & 1 == 1 {
            result *= base_u64;
            result %= mod_u64;
        }
        exp_u64 >>= 1;
        base_u64 *= base_u64;
        base_u64 %= mod_u64;
    }

    return T::from_u64(result).unwrap();
    //return result as u128;
}

pub fn is_invertible(f: &Vec<i32>, p: u128) -> bool {
    // asserts if the polynomial f is invertible mod X^n + 1
    // case for p=2 works because in integers mod 2, a polynomial is invertible <->
    // sum of coefficients is odd <-> non-zero constant-term
    if p == 2{
        let mut sum: i32 = 0;
        for i in 0..f.len() {
            sum += f[i];
            sum %= 2;
        }
        return sum == 1;
    }
    // if p is some other prime, we can use NTT representation of f to check invertibility

    return false
}

pub fn l2norm(f: &Vec<i32>) -> i64 {
    // returns the l2 norm of polynomial/vector f as f[0]^2 + f[1]^2 +..+ f[n]^2
    let mut sum: i64 = 0;
    for i in 0..f.len() {
        sum += (f[i] as i64).pow(2);
    }
    return sum;
}

pub fn adjoint(f: &Vec<i32>) -> Vec<i32> {
    // computes the (hermitian) adjoint of a polynomial f 
    let mut fstar = f.clone();
    for i in 1..f.len() {
        fstar[i] = -f[f.len() - i];
    }
    return fstar;
}

pub fn poly_add(f: &Vec<i32>, g: &Vec<i32>) -> Vec<i32> {
    // performs standard polynomial addition of two polynomials with mod p
    assert_eq!(f.len(), g.len());
    let mut q = vec![0; f.len()];
    for i in 0..q.len() {
        q[i] = f[i] + g[i];
    }
    return q;
}

pub fn poly_mult(f: &Vec<i32>, g: &Vec<i32>, p:i32) -> Vec<i32> {
    // performs standard polynomial multiplication of two polynomials with mod p
    let mut q = vec![0; f.len() + g.len() - 1];

    for i in 0..f.len() {
        for j in 0..g.len() {
            q[i + j] += f[i] * g[j];
            q[i+j] %= p;
        }
    }
    return q;
}
