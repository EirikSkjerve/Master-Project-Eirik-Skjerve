use crate::hawkverify::ptc;
use crate::ntt::*;
use crate::ntt_constants::brv;
use crate::utils::{adjoint, mod_pow, modulo};
use num_complex::Complex;
use std::f64::consts::PI;

pub fn delta(k: usize) -> (i64, i64) {
    let i = Complex::new(0.0, 1.0);
    let pi_complex = Complex::new(PI, 0.0);
    let temp = Complex::new(brv(k as u128, 10) as f64, 0.0);
    let const_complex = Complex::new(2048.0, 0.0);

    let temp2: Complex<f64> = (2.0 * i * pi_complex * temp) / const_complex;

    let mut d = temp2.exp();

    let base: f64 = 2.0;
    let factor: Complex<f64> = Complex::new(base.powi(31), 0.0);
    d *= factor;

    let re: i64 = d.re.round() as i64;
    let im: i64 = d.im.round() as i64;

    (re, im)
}

fn sign(x: i64) -> i64 {
    // returns the sign of an i64:
    // 1 if x is negative, 0 otherwise
    if x < 0 {
        return 1;
    }
    0
}

pub fn rebuildw0(
    q00: &Vec<i64>,
    q01: &Vec<i64>,
    w1: &Vec<i64>,
    h0: &Vec<i64>,
    highs0: usize,
    highs1: usize,
    high00: usize,
    high01: usize,
) -> Vec<i64> {
    let n = q00.len();

    let base_i64: i64 = 2;

    let cw1 = 1 << (29 - (1 + highs1));
    let cq00 = 1 << (29 - high00);
    let cq01 = 1 << (29 - high01);

    let cs0 = (2 * (cw1) * (cq01)) / (n as i64 * cq00);
    let w1_fft = fft(&ptc(w1, cw1));

    let mut z00 = q00.clone();

    // some test here
    if z00[0] < 0 {
        // TODO return None here instead
        println!("rebuild failure 1");
        return vec![0];
    }

    z00[0] = 0;

    // compute fft(q00*cq00) and fft(q01*cq01)
    let q00_fft = fft(&ptc(&z00, cq00));
    let mut q01_fft = fft(&ptc(&q01, cq01));

    let alpha = (2 * (cq00) * (q00[0])) / n as i64;

    let n_uz = n as usize;
    for u in 0..(n_uz / 2) {
        let mut x_re = q01_fft[u] * w1_fft[u];
        x_re -= q01_fft[u + (n_uz / 2)] * w1_fft[u + (n_uz / 2)];

        let mut x_im = q01_fft[u] * w1_fft[u + (n_uz / 2)];
        x_im += q01_fft[u + (n_uz / 2)] * w1_fft[u];

        let (x_re, z_re) = (x_re.abs(), sign(x_re));
        let (x_im, z_im) = (x_im.abs(), sign(x_im));

        let v = alpha + q00_fft[u];

        if v <= 0
            || v >= base_i64.pow(32)
            || (x_re) >= v * base_i64.pow(32)
            || (x_im) >= v * base_i64.pow(32)
        {
            // return None here instead
            return vec![0];
        }

        // let y_re = x_re.div_floor(&v);
        let y_re = num::Integer::div_floor(&x_re, &v);
        // let y_im = x_im.div_floor(&v);
        let y_im = num::Integer::div_floor(&x_im, &v);

        q01_fft[u] = y_re - (2 * z_re * y_re);
        q01_fft[u + (n_uz / 2)] = y_im - (2 * z_im * y_im);
    }

    let t = ifft(&q01_fft);

    let mut w0: Vec<i64> = vec![0; n_uz];

    for u in 0..n_uz {
        let v = (cs0 * h0[u]) + t[u];
        // let z = (v + cs0).div_floor(&(2 * cs0));
        let z = num::Integer::div_floor(&(v + cs0), &(2 * cs0));

        if z < -base_i64.pow(highs0 as u32) || z >= base_i64.pow(highs0 as u32) {
            // return None here instead
            println!("rebuild failure 3");
            return vec![0];
        }

        w0[u] = h0[u] - (2 * z);
    }

    w0
}

pub fn fft(f: &Vec<i64>) -> Vec<i64> {
    let n = f.len();
    let mut f_fft: Vec<i64> = f.clone();
    let mut t = n / 2;
    let mut m = 2;

    while m < n {
        let mut v0 = 0;
        for u in 0..(m / 2) {
            let e = delta(u + m);
            let e_re: i64 = e.0 as i64;
            let e_im: i64 = e.1 as i64;

            for v in v0..v0 + (t / 2) {
                let x1_re: i64 = f_fft[v];
                let x1_im: i64 = f_fft[v + (n / 2)];

                let x2_re: i64 = f_fft[v + (t / 2)];
                let x2_im: i64 = f_fft[v + (t / 2) + (n / 2)];

                let t_re = (x2_re * e_re) - (x2_im * e_im);
                let t_im = (x2_re * e_im) + (x2_im * e_re);

                f_fft[v] = ((x1_re << 31) + t_re) >> 32;
                f_fft[v + (n / 2)] = ((x1_im << 31) + t_im) >> 32;
                f_fft[v + (t / 2)] = ((x1_re << 31) - t_re) >> 32;
                f_fft[v + (t / 2) + (n / 2)] = ((x1_im << 31) - t_im) >> 32;
            }
            v0 += t;
        }
        t /= 2;
        m *= 2;
    }

    f_fft
}

pub fn ifft(f_fft: &Vec<i64>) -> Vec<i64> {
    let n = f_fft.len();

    let mut f = f_fft.clone();
    let mut t = 2;
    let mut m = n / 2;

    while m > 1 {
        let mut v0 = 0;
        for u in 0..(m / 2) {
            let e = delta(u + m);
            let e_re = e.0;
            let e_im = -e.1;

            for v in v0..v0 + (t / 2) {
                let x1_re = f[v];
                let x1_im = f[v + (n / 2)];

                let x2_re = f[v + (t / 2)];
                let x2_im = f[v + (t / 2) + (n / 2)];

                let t1_re = x1_re + x2_re;
                let t1_im = x1_im + x2_im;

                let t2_re = x1_re - x2_re;
                let t2_im = x1_im - x2_im;

                f[v] = t1_re >> 1;
                f[v + (n / 2)] = t1_im >> 1;
                f[v + (t / 2)] = (t2_re * e_re - t2_im * e_im) >> 32;
                f[v + (t / 2) + (n / 2)] = (t2_re * e_im + t2_im * e_re) >> 32;
            }
            v0 += t;
        }
        t *= 2;
        m /= 2;
    }

    f
}

pub fn poly_qnorm(q00: &Vec<i64>, q01: &Vec<i64>, w0: &Vec<i64>, w1: &Vec<i64>, p: i64) -> i64 {
    // compute q-norm of w = (w0, w1)

    let n = q00.len();

    // compute ntt representation of q and w
    let q00ntt = ntt(q00.clone(), p);
    let q01ntt = ntt(q01.clone(), p);
    let w0ntt = ntt(w0.clone(), p);
    let w1ntt = ntt(w1.clone(), p);

    let mut d: Vec<i64> = Vec::with_capacity(n);
    for i in 0..n {
        let val = modulo(w1ntt[i] * mod_pow(q00ntt[i], p - 2, p), p);
        d.push(val);
    }

    let w1_adj_ntt = ntt(adjoint(&w1), p);
    let mut c: Vec<i64> = Vec::with_capacity(n);
    for i in 0..n {
        let val = modulo(d[i] * w1_adj_ntt[i], p);
        c.push(val);
    }

    let mut acc: i64 = c.iter().sum();

    let mut e: Vec<i64> = Vec::with_capacity(n);
    for i in 0..n {
        let val = modulo(w0ntt[i] + (d[i] * q01ntt[i]), p);
        e.push(val);
    }

    let e_adj_ntt = nttadj(&e, p);
    for i in 0..n {
        // preprocessing
        let temp1 = modulo(q00ntt[i], p) as i64;
        let temp2 = modulo(e[i], p) as i64;
        let temp3 = modulo(e_adj_ntt[i], p) as i64;
        let val = modulo(temp1 as i128 * temp2 as i128 * temp3 as i128, p as i128) as i64;

        acc = modulo(acc, p);
        acc += val;
    }

    modulo(acc, p)
}
