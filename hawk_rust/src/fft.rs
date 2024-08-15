use std::f64::consts::PI;

use libm::Libm;
use num::Zero;
use num_complex::Complex;
use crate::fft_constants;

pub fn fft(f: &Vec<f64>) -> Vec<Complex<f64>> {
    let n = f.len();

    // base case
    if n == 2 {
        let mut f_fft: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); n];
        f_fft[0] = Complex::new(f[0], f[1]);
        f_fft[1] = Complex::new(f[0], -f[1]);
        return f_fft;
    }

    // getting roots

    let w = fft_constants::get_roots(n);

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

pub fn ifft(f_fft: &Vec<Complex<f64>>) -> Vec<f64> {
    let n = f_fft.len();

    // base case
    if n == 2 {
        let mut f: Vec<f64> = vec![0.0; n];
        f[0] = f_fft[0].re;
        f[1] = f_fft[0].im;
        return f;

    }

    let mut f: Vec<f64> = vec![0.0; n];

    let w = fft_constants::get_roots(n);
    let m = n/2;

    let mut f0_fft: Vec<Complex<f64>> = vec![Complex::new(0.0,0.0); m];
    let mut f1_fft: Vec<Complex<f64>> = vec![Complex::new(0.0,0.0); m];

    for i in 0..m {
        f0_fft[i] = 0.5 * (f_fft[2*i] + f_fft[2*i+1]);
        f1_fft[i] = 0.5 * (f_fft[2*i] - f_fft[2*i+1])* w[2*i].conj();

    }

    // recursively call ifft on halves
    let f0 = ifft(&f0_fft);
    let f1 = ifft(&f1_fft);

    // merge the splits
    for i in 0..m {
        f[2*i] = f0[i];
        f[2*i+1] = f1[i];
    }


    return f;
}

// fft where input is i64: converts the polynomial into complex domain
pub fn f_fft(f: &Vec<i64>) -> Vec<Complex<f64>> {
    let mut f_f: Vec<f64> = Vec::new();
    for i in 0..f.len() {
        f_f.push(f[i] as f64);
    }

    return fft(&f_f);
}
