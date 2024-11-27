use crate::ntt::{intt, ntt};
use num::traits::{FromPrimitive, Num, PrimInt};

use num_bigint::{BigInt, ToBigInt};
use num_traits::{Signed, ToPrimitive};

pub fn bytes_to_poly(h: &[u8], n: usize) -> Vec<i64> {
    /*
     * converts a byte-array to a vector/polynomial
     */
    let mut res: Vec<i64> = vec![0; n];
    for i in 0..n {
        res[i] = ((h[i / 8] >> (i % 8)) & 1) as i64;
    }
    return res;
}

// implements integer modulation for positive and negative integers
pub fn modulo<T: PrimInt>(a: T, b: T) -> T
where
    T: Num + FromPrimitive + std::fmt::Display,
{
    //
    // returns a mod b
    //
    // example: modulo(15, 4) -> 3
    //          modulo(-5, 3) -> 1
    //
    // convert the inputs to u128
    //

    // if a > T::from(0).unwrap() {
    //     let a_u128 = a.to_u128().unwrap();
    //     let b_u128 = b.to_u128().unwrap();
    //
    //     // perform the calculations
    //     let result = ((a_u128 % b_u128) + b_u128) % b_u128;
    //
    //     return T::from_u128(result).unwrap();
    // }

    let a_i128 = a.to_i128().unwrap();
    let b_i128 = b.to_i128().unwrap();

    // perform the calculations
    let result = ((a_i128 % b_i128) + b_i128) % b_i128;

    T::from_i128(result).unwrap()
}

pub fn bigint_vec(v: &Vec<i64>) -> Vec<BigInt> {
    //
    // returns a copy of input vector with each entry converted to type BigInt.
    // Signs are preserved
    //

    let mut v_big: Vec<BigInt> = Vec::new();
    for i in v.iter() {
        v_big.push(i.to_bigint().unwrap());
    }

    v_big
}

pub fn bigint_to_f64_vec(a: Vec<BigInt>) -> Vec<f64> {
    //
    // returns a copy of input BigInt vector with each entry converted to type f64
    //

    let n = a.len();
    let mut res: Vec<f64> = Vec::with_capacity(n);

    for i in 0..n {
        if let Some(res_i) = &a[i].to_f64() {
            res.push(*res_i);
        } else {
            println!("Could not convert to float");
        }
    }

    res
}

pub fn bigint_to_i64_vec(a: Vec<BigInt>) -> Vec<i64> {
    //
    // returns a copy of input BigInt vector with each entry converted to type i64
    //

    let n = a.len();
    let mut res: Vec<i64> = Vec::with_capacity(n);

    for i in 0..n {
        if let Some(res_i) = &a[i].to_i64() {
            res.push(*res_i);
        } else {
            println!("Could not convert to i64");
        }
    }

    res
}
// implements fast binary exponentiation for computing base^exp mod modulus
// inputs base, exponent and modulus as generic, and returns a u128
// note that this is an already implemented method for BigInt, which is probably better to use
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

    T::from_u128(result).unwrap()
}

pub fn is_invertible(f: &Vec<i64>, p: i64) -> bool {
    //
    // asserts if the polynomial f is invertible mod X^n + 1
    // case for p=2 works because in integers mod 2, a polynomial is invertible <->
    // sum of coefficients is odd <-> non-zero constant-term
    //

    if p == 2 {
        let mut sum: i64 = 0;
        for i in 0..f.len() {
            sum += f[i];
            sum = modulo(sum, 2);
        }
        return sum == 1;
    }

    // if p is some other prime, we use NTT representation of f to check invertibility
    let f_ntt = ntt(f.clone(), p);
    for i in 0..f.len() {
        if f_ntt[i] == 0 {
            return false;
        }
    }

    true
}

pub fn l2norm(f: &Vec<i64>) -> i64 {
    //
    // returns the l2 norm of polynomial/vector f as f[0]^2 + f[1]^2 +..+ f[n]^2
    //

    let mut sum: i64 = 0;
    for i in 0..f.len() {
        sum += f[i] * f[i];
    }
    sum
}

pub fn infnorm(f: &Vec<i64>) -> i64 {
    //
    // returns the inf norm of polynomial/vector f as max(f[0], f[1], ..., f[n])
    //

    let max = f.iter().map(|x| x.abs()).max().unwrap();
    max
}

pub fn adjoint(f: &Vec<i64>) -> Vec<i64> {
    //
    // computes the hermitian adjoint of a polynomial f
    //

    let mut fstar = f.clone();
    for i in 1..f.len() {
        fstar[i] = -f[f.len() - i];
    }
    fstar
}

pub fn poly_add(f: &Vec<i64>, g: &Vec<i64>) -> Vec<i64> {
    //
    // performs standard polynomial addition of two polynomials
    //

    assert_eq!(f.len(), g.len());
    let mut q = vec![0; f.len()];
    for i in 0..q.len() {
        q[i] = f[i] + g[i];
    }
    q
}

pub fn poly_sub(f: &Vec<i64>, g: &Vec<i64>) -> Vec<i64> {
    //
    // performs standard polynomial subtraction of two polynomials
    //

    assert_eq!(f.len(), g.len());
    let mut q = vec![0; f.len()];
    for i in 0..q.len() {
        q[i] = f[i] - g[i];
    }
    q
}

pub fn poly_mult_ntt(f: &Vec<i64>, g: &Vec<i64>, p: i64) -> Vec<i64> {
    //
    // perform polynomial multiplication of f and g mod p mod X^n + 1
    // using the Number Theoretical Transform for efficiency
    //

    // length of f and g should be the same
    assert_eq!(f.len(), g.len());

    let n = f.len();

    // ntt representation of f and g
    let f_ntt = ntt(f.clone(), p);
    let g_ntt = ntt(g.clone(), p);

    // polynomials in ntt representation can be multiplied point-wise
    let mut fg_ntt: Vec<i64> = vec![0; n];
    for i in 0..n {
        fg_ntt[i] = modulo(f_ntt[i] * g_ntt[i], p as i64);
    }

    // convert f*g back to normal polynomial form
    let mut fg = intt(fg_ntt, p);

    // center the polynomial around 0 / adjust coefficients
    for i in 0..n {
        if fg[i] > (p - 1) / 2 {
            fg[i] -= p;
        }
    }

    fg
}

fn transpose(f: &Vec<Vec<i64>>) -> Vec<Vec<i64>> {
    // transposes a 2d list, i.e. swaps every i and j position

    let mut ft: Vec<Vec<i64>> = Vec::new();

    for i in 0..f.len() {
        let mut temp: Vec<i64> = Vec::new();
        for j in 0..f[i].len() {
            temp.push(f[j][i]);
        }
        ft.push(temp);
    }

    ft
}

pub fn rot(f: &Vec<i64>) -> Vec<Vec<i64>> {
    // implements the rot() functionality from Hawk spec. paper
    // rot(f) = (vec(f), vec(fX),..., vec(fX^{n-1}))

    let n = f.len();
    let mut rotf: Vec<Vec<i64>> = Vec::new();
    rotf.push(f.clone());
    let mut cur = f.clone();
    for i in 1..n {
        let temp = -cur.remove(n - 1);
        cur.insert(0, temp);
        rotf.push(cur.clone());
    }

    rotf
}

pub fn rot_key(f: &Vec<i64>, g: &Vec<i64>, bigf: &Vec<i64>, bigg: &Vec<i64>) -> Vec<Vec<i64>> {
    // performs rot on a key, either B or Q
    // so for example rot_key(B) = (rot(f), rot(F), rot(g), rot(G)), and likewise for q
    // assumes all input polynomials are of the same length

    let n = f.len();

    // rotate each part individually
    let rotf = rot(f);
    let rotg = rot(g);
    let rotbigf = rot(bigf);
    let rotbigg = rot(bigg);

    // result vector
    let mut res: Vec<Vec<i64>> = vec![Vec::with_capacity(2 * n); 2 * n];

    // join rot(f) and rot(F) and add the resulting rows to res
    for i in 0..n {
        let mut rowi: Vec<i64> = Vec::with_capacity(2 * n);
        for j in 0..n {
            rowi.push(rotf[j][i]);
        }
        for j in 0..n {
            rowi.push(rotbigf[j][i]);
        }
        res[i] = rowi;
    }

    // join rot(g) and rot(G) and add the resulting rows to res
    for i in n..2 * n {
        let mut rowi: Vec<i64> = Vec::with_capacity(2 * n);
        for j in 0..n {
            rowi.push(rotg[j][i - n]);
        }
        for j in 0..n {
            rowi.push(rotbigg[j][i - n]);
        }
        res[i] = rowi;
    }

    res
}
