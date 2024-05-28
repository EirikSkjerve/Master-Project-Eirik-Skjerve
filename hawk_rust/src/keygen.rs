use crate::rngcontext::{shake256x4, RngContext};
use crate::utils::{is_invertible, l2norm, adjoint};

pub fn hawkkeygen(logn: u16, rng: Option<RngContext>) {
    // checks if rng-context is initialized or not. If not, initialize a new one and recusively call hawkkeygen
    let mut rng = match rng {
        Some(rng) => rng,
        None => {
            let new_rng = RngContext::new(1337); // this should be an actual random number
            return hawkkeygen(logn, Some(new_rng));
        }
    };

    // generate f and g
    let f_g = generate_f_g(rng.rnd(128) as usize, logn);
    let f = f_g.0.clone();
    let g = f_g.1.clone();

    // checks if f and g are invertible mod X^n + 1 and mod 2
    // if not, restart
    if !(is_invertible(&f) && is_invertible(&g)) {
        return hawkkeygen(logn, Some(rng));
    }

    let n = 1 << logn;

    // checks if the norm of f and g is large enough
    // if not, restart
    if ((l2norm(&f) + l2norm(&g)) as f64) <= 2.0 * (n as f64) * (1.042) {
        return hawkkeygen(logn, Some(rng));
    }

    // construct the adjoints of f and g
    let fstar = adjoint(&f);
    let gstar = adjoint(&g);

    println!("f: {:?}, \n g: {:?}", f, g);

}

// generates polynomials f and g
fn generate_f_g(seed: usize, logn: u16) -> (Vec<i32>, Vec<i32>) {
    // expand logn -> n
    let n = 1 << logn;
    let b = n / 64;
    assert!(b == 4 || b == 8 || b == 16);

    // get an array of 64-bit values
    let y = shake256x4(&seed.to_ne_bytes(), 2 * n * b / 64);

    // construct a sequence of bits from y
    let mut ybits: Vec<u8> = vec![0; b * 2 * n];
    for (j, y) in y.iter().enumerate() {
        for bi in 0..64 {
            ybits[j * 64 + bi] = ((y >> bi) & 1) as u8;
        }
    }

    // generate f and g from centered binomial distribution
    // if e.g. n = 256, b = 4, so f and g consists of random numbes from the interval [-2,-1,0,1,2]
    // if n = 512, b = 8, inteval = [-4,...,4]
    // if n = 1024, b = 16, interval = [-8,...,8]
    let mut f: Vec<i32> = vec![0; n];
    let mut sum: i32 = 0;

    // get some bounded number from ybits
    for i in 0..n {
        sum = 0;
        for j in 0..b {
            sum += ybits[i * b + j] as i32;
        }

        // center the number around 0
        f[i] = sum - (b / 2) as i32;
    }

    // reuse "sum" variable here
    let mut g: Vec<i32> = vec![0; n];

    for i in 0..n {
        sum = 0;
        for j in 0..b {
            sum += ybits[i + n * b + j] as i32;
        }

        g[i] = sum - (b / 2) as i32;
    }

    // returns the two vectors
    return (f, g);
}
