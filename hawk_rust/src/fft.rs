use std::f64::consts::PI;

use num::Zero;
use num_complex::Complex;
use libm::Libm;


fn fft(f: &Vec<Complex<f64>>) -> Vec<Complex<f64>> {
    let n = f.len();

    // base case 
    if n == 2 {
        // let mut y: Vec<Complex<f64>> = vec![Complex::new(0.0,0.0); n];
        return f.to_vec();
    }

    // calculating roots
    let theta = -2.0*PI/(n as f64);
    let mut w: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); n];

    let mut temp: f64 = 0.0;
    for i in (0..n) {
        temp = theta* i as f64;
        w[i] = Complex::new(Libm::<f64>::cos(temp), Libm::<f64>::sin(temp));
    }

    // initialize empty odd/even vectors
    let mut f_even: Vec<Complex<f64>> = Vec::new(); 
    let mut f_odd: Vec<Complex<f64>> = Vec::new();

    // split f into two
    for i in 0..n {
        if i%2==0{
            f_even.push(f[i]);
        }
        else{
            f_odd.push(f[i]);
        }
    }

    // recursively call fft on halves
    let y_even = fft(&f_even);
    let y_odd = fft(&f_odd);

    // initialize empty result vector
    let mut y: Vec<Complex<f64>> = vec![Complex::new(0.0,0.0); n];

    let m = n/2;

    let mut temp = Complex::zero();
    for i in 0..m {

        temp = w[2*i] * y_odd[i];
        y[2*i] = y_even[i] + temp;
        y[2*i+1] = y_even[i] - temp;   
    }
    /*
    // initialize working variables
    let mut w_yodd_k = Complex::zero();
    let mut yeven_k = Complex::zero();

    // perform the calculations
    for k in 0..m {

       w_yodd_k = w[k] * y_odd[k];
       yeven_k = y_even[k];

       y[k] = yeven_k + w_yodd_k;
       y[k + m] = yeven_k - w_yodd_k;
    }
    */

    // return the result
    return y;

}

// fft where input is i64: converts the polynomial into complex domain
pub fn f_fft(f: &Vec<i64>) -> Vec<Complex<f64>> {
    
    let poly = f.clone();

    let poly_complex = poly.into_iter()
        .map(|x| Complex::new(x as f64, 0.0))
        .collect();

    return fft(&poly_complex);

}
