use crate::fft_constants;
use num_complex::{Complex, ComplexFloat};

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

    // join polynomials
    for i in 0..m {
        let temp: Complex<f64> = w[2 * i] * y_o[i];
        y[2 * i] = y_e[i] + temp;
        y[2 * i + 1] = y_e[i] - temp;
    }

    // return the result
    return y;
}

pub fn fft_i64(f: &Vec<i64>) -> Vec<Complex<f64>> {
    let f_f: Vec<f64> = f.iter().map(|&x| x as f64).collect();
    return fft(&f_f);
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
    let m = n / 2;

    let mut f0_fft: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); m];
    let mut f1_fft: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); m];

    for i in 0..m {
        f0_fft[i] = 0.5 * (f_fft[2 * i] + f_fft[2 * i + 1]);
        f1_fft[i] = 0.5 * (f_fft[2 * i] - f_fft[2 * i + 1]) * w[2 * i].conj();
    }

    // recursively call ifft on halves
    let f0 = ifft(&f0_fft);
    let f1 = ifft(&f1_fft);

    // merge the splits
    for i in 0..m {
        f[2 * i] = f0[i];
        f[2 * i + 1] = f1[i];
    }

    return f;
}

pub fn inverse_fft(p: &Vec<i64>) -> Vec<f64> {
    // convert to vector of floats/rationals
    let mut p_f: Vec<f64> = vec![0.0; p.len()];
    for i in 0..p.len() {
        p_f[i] = p[i] as f64;
    }

    let m = p.len() / 2;
    let mut p_fft = fft(&p_f);
    //println!("p_fft: {:?}", p_fft);

    for u in 0..m {
        p_fft[u].re = 1.0 / p_fft[u].re;
        p_fft[u + m].re = 0.0;

        // doing this to match the reference implementation
        p_fft[u].im = 0.0;
        p_fft[u + m].im = 0.0;
    }

    return ifft(&p_fft);
}

pub fn mul_fft(f: &Vec<Complex<f64>>, g: &Vec<Complex<f64>>) -> Vec<Complex<f64>> {
    let n = f.len();
    let mut res: Vec<Complex<f64>> = Vec::with_capacity(n);

    for i in 0..n {
        res.push(f[i] * g[i]);
    }

    return res;
}

pub fn div_fft(f: &Vec<Complex<f64>>, g: &Vec<Complex<f64>>) -> Vec<Complex<f64>> {
    let n = f.len();
    let mut res: Vec<Complex<f64>> = Vec::with_capacity(n);

    for i in 0..n {
        res.push(f[i] / g[i]);
    }

    return res;
}

pub fn add_fft(f: &Vec<Complex<f64>>, g: &Vec<Complex<f64>>) -> Vec<Complex<f64>> {
    let n = f.len();
    let mut res: Vec<Complex<f64>> = Vec::with_capacity(n);

    for i in 0..n {
        res.push(f[i] + g[i]);
    }

    return res;
}

pub fn adj_fft(f: &Vec<Complex<f64>>) -> Vec<Complex<f64>> {
    let n = f.len();
    let mut res: Vec<Complex<f64>> = Vec::with_capacity(n);

    for i in 0..n {
        res.push(f[i].conj());
    }

    return res;
}
