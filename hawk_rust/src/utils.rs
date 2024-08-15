use crate::ntt_fft::{intt, ntt};
use num::traits::{FromPrimitive, Num, PrimInt};

pub fn bin(a: u128, x: usize) -> Vec<u8> {
    /*
    Converts an integer to binary representation in a vector of arbitrary size
    Call bin(a, x=0) for default binary size
    */
    let mut res: Vec<u8> = vec![];
    let mut b = a;

    if x == 0 {
        // returns a as binary with *default* size
        return bin(a, (a as f64).log2().ceil() as usize + 1);
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

// implements integer modulation
pub fn modulo<T: PrimInt>(a: T, b: T) -> T
where
    T: Num + FromPrimitive,
{
    // convert the inputs to u64
    let a_i64 = a.to_i64().unwrap();
    let b_i64 = b.to_i64().unwrap();

    // perform the calculations
    let result = ((a_i64 % b_i64) + b_i64) % b_i64;

    return T::from_i64(result).unwrap();
}

// implements fast binary exponentiation for computing base^exp mod modulus
// inputs base, exponent and modulus as generic, and returns a u128
pub fn mod_pow<T: PrimInt>(base: T, exp: T, modulus: T) -> T
where
    T: Num + FromPrimitive,
{
    // convert the inputs to u64
    let mut base_u128 = base.to_u128().unwrap();
    let mut exp_u128 = exp.to_u128().unwrap();
    let mod_u128 = modulus.to_u128().unwrap();

    // perform the algorithm
    let mut result = 1;
    base_u128 %= mod_u128;
    while exp_u128 > 0 {
        if exp_u128 & 1 == 1 {
            result *= base_u128;
            result %= mod_u128;
        }
        exp_u128 >>= 1;
        base_u128 *= base_u128;
        base_u128 %= mod_u128;
    }

    return T::from_u128(result).unwrap();
}

pub fn is_invertible(f: &Vec<i64>, p: u32) -> bool {
    // asserts if the polynomial f is invertible mod X^n + 1
    // case for p=2 works because in integers mod 2, a polynomial is invertible <->
    // sum of coefficients is odd <-> non-zero constant-term
    if p == 2 {
        let mut sum: i64 = 0;
        for i in 0..f.len() {
            sum += f[i];
            sum = modulo(sum, 2);
        }
        return sum == 1;
    }
    // if p is some other prime, we can use NTT representation of f to check invertibility
    let f_ntt = ntt(f.clone(), p);
    for i in 0..f.len() {
        if f_ntt[i] == 0 {
            return false;
        }
    }

    return true;
}

pub fn l2norm(f: &Vec<i64>) -> i64 {
    // returns the l2 norm of polynomial/vector f as f[0]^2 + f[1]^2 +..+ f[n]^2
    let mut sum: i64 = 0;
    for i in 0..f.len() {
        sum += (f[i]).pow(2);
    }
    println!("{}", sum);
    return sum;
}

pub fn adjoint(f: &Vec<i64>) -> Vec<i64> {
    // computes the (hermitian) adjoint of a polynomial f
    let mut fstar = f.clone();
    for i in 1..f.len() {
        fstar[i] = -f[f.len() - i];
    }
    return fstar;
}

pub fn poly_add(f: &Vec<i64>, g: &Vec<i64>) -> Vec<i64> {
    // performs standard polynomial addition of two polynomials with mod p
    assert_eq!(f.len(), g.len());
    let mut q = vec![0; f.len()];
    for i in 0..q.len() {
        q[i] = f[i] + g[i];
    }
    return q;
}

pub fn poly_mult(f: &Vec<i64>, g: &Vec<i64>, p: i64) -> Vec<i64> {
    // performs standard polynomial multiplication of two polynomials with mod p
    let mut q = vec![0; f.len() + g.len() - 1];

    for i in 0..f.len() {
        for j in 0..g.len() {
            q[i + j] += f[i] * g[j];
            q[i + j] %= p;
        }
    }
    return q;
}

pub fn poly_mult_ntt(f: &Vec<i64>, g: &Vec<i64>, p: u32) -> Vec<i64> {
    // length of f and g should be the same
    assert_eq!(f.len(), g.len());

    let n = f.len();

    // ntt representation of f and g
    let f_ntt = ntt(f.clone(), p);
    let g_ntt = ntt(g.clone(), p);

    let mut fg_ntt: Vec<i64> = vec![0; n];
    for i in 0..n {
        fg_ntt[i] = modulo(f_ntt[i] * g_ntt[i], p as i64);
    }

    let mut fg = intt(fg_ntt, p);

    for i in 0..n {
        if fg[i] > (p as i64 - 1) / 2 {
            fg[i] -= p as i64;
        }
    }

    return fg;
}
