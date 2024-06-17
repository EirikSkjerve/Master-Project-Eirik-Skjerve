pub fn bin(a: u128, x: usize) -> Vec<u8> {
    /*
    Converts an integer to binary representation in a vector of arbitrary size
    */
    let mut res: Vec<u8> = vec![];
    let mut b = a;

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

pub fn is_invertible(f: &Vec<i32>, p: u128) -> bool {
    // asserts if the polynomial f is invertible mod X^n + 1
    if p == 2{
    let mut sum: i32 = 0;
    for i in 0..f.len() {
        sum += f[i];
        sum %= 2;
    }
    return sum == 1;
    }

    return false
}

pub fn l2norm(f: &Vec<i32>) -> i64 {
    let mut sum: i64 = 0;
    for i in 0..f.len() {
        sum += (f[i] as i64).pow(2);
    }
    return sum;
}

pub fn adjoint(f: &Vec<i32>) -> Vec<i32> {
    let mut fstar = f.clone();
    for i in 1..f.len() {
        fstar[i] = -f[f.len() - i];
    }
    return fstar;
}

// performs polynomial addition of two polynomials without modulation
pub fn poly_add(f: &Vec<i32>, g: &Vec<i32>) -> Vec<i32> {
    assert_eq!(f.len(), g.len());
    let mut q = vec![0; f.len()];
    for i in 0..q.len() {
        q[i] = f[i] + g[i];
    }
    return q;
}

// performs polynomial multiplication of two polynomials without modulation
pub fn poly_mult(f: &Vec<i32>, g: &Vec<i32>) -> Vec<i32> {
    let mut q = vec![0; f.len() + g.len() - 1];

    for i in 0..f.len() {
        for j in 0..g.len() {
            q[i + j] += f[i] * g[j];
        }
    }
    return q;
}