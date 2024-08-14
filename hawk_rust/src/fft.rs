use std::f64::consts::PI;

use libm::Libm;
use num::Zero;
use num_complex::Complex;
use crate::fft_constants;

fn fft(f: &Vec<f64>) -> Vec<Complex<f64>> {
    println!("f: {:?}", f);
    let n = f.len();

    // base case
    if n == 2 {
        let mut f_fft: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); n];
        f_fft[0] = Complex::new(f[0], f[1]);
        f_fft[1] = Complex::new(f[0], -f[1]);
        return f_fft;
    }

    // getting roots

    println!("w: {:?}", w);

    let mut f_e: Vec<f64> = Vec::new();
    let mut f_o: Vec<f64> = Vec::new();

    // split f into two
    for i in 0..n {
        if i % 2 == 0 {
            f_e.push(f[i]);
        } else {
            f_o.push(f[i]);
        }
    }

    // recursively call fft on halves
    let y_e = fft(&f_e);
    let y_o = fft(&f_o);

    // initialize empty result vector
    let mut y: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); n];

    let m = n / 2;

    let mut temp = Complex::zero();

    // join polynomials
    for i in 0..m {
        temp = w[2 * i] * y_o[i];
        y[2 * i] = y_e[i] + temp;
        y[2 * i + 1] = y_e[i] - temp;
    }

    // return the result
    return y;
}

// fft where input is i64: converts the polynomial into complex domain
pub fn f_fft(f: &Vec<i64>) -> Vec<Complex<f64>> {
    let mut f_f: Vec<f64> = Vec::new();
    for i in 0..f.len() {
        f_f.push(f[i] as f64);
    }

    return fft(&f_f);
}
