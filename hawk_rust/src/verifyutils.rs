use num_complex::Complex;
use std::f64::consts::PI;
use crate::ntt::*;
use crate::utils::{adjoint, mod_pow, modulo};


pub fn delta(k: usize) -> (u32, u32) {

    let i = Complex::new(0.0, 1.0);
    let pi_complex = Complex::new(PI, 0.0);
    let temp = Complex::new(brv(k as u32, 10) as f64, 0.0);
    let const_complex = Complex::new(2048.0, 0.0);

    let temp2: Complex<f64> = (2.0*i * pi_complex * temp) / const_complex; 

    let mut d = temp2.exp();

    let base: f64 = 2.0;
    let factor: Complex<f64> = Complex::new(base.powi(31), 0.0);
    d *= factor;

    let re: u32 = d.re.round() as u32;
    let im: u32 = d.im.round() as u32;

    return (re, im)

}

// implementation of bit reversal of an integer
fn brv(x: u32, log_n: usize) -> u32 {
    // following is HAWK-teams implementation. Might go about it a different way
    // no idea how it works though

    let mut r = 0;
    for i in 0..log_n {
        r |= ((x >> i) & 1) << (log_n - 1 - i);
    }
    return r;
}


fn sign(x: i64) -> i64 {
    if x < 0{
        return 1;
    }
    return 0;
}

fn scale_vec(f: &Vec<i32>, c: i32) -> Vec<i32> {
    let mut f_scaled = Vec::with_capacity(f.len());
    for i in 0..f.len() {
        println!("{} * {} =...", f[i], c);
        f_scaled.push(f[i]*c);
    }
    return f_scaled;
}

pub fn rebuildw0(logn: usize, q00: &Vec<i32>, q01: &Vec<i32>, w1: &Vec<i32>, h0: &Vec<i32>) -> Vec<i32> {

    let n = 1 << logn;
    let highs0 = 12;
    let highs1 = 9;
    let high00 = 9;
    let high01 = 11;

    let base: i32 = 2;
    let base_i64: i64 = 2; 
    let cw1 = base.pow(29 - (1 + highs1));
    let cq00 = base.pow(29 - high00);
    let cq01 = base.pow(29 - high01);
    let cs0 = ((2*(cw1 as i64)*(cq01 as i64)) / ((n*cq00) as i64)) as i32;

    println!("here 1");
    let w1_fft = fft(&scale_vec(&w1, cw1));

    let mut z00 = q00.clone();
    
    // some test here
    if z00[0] < 0 {
        return vec![0];
    }

    z00[0] = 0;

    println!("here 2");
    let q00_fft = fft(&scale_vec(&z00, cq00));
    println!("here 3");
    let scaled_q01 = scale_vec(&q01, cq01);
    println!("here 4");
    let mut q01_fft = fft(&scale_vec(&q01, cq01));

    let alpha = (2*cq00*q00[0])/n;

    let n_uz = n as usize;
    for u in 0..(n_uz/2) {
        let mut x_re = (q01_fft[u] * w1_fft[u+(n_uz/2)]) as i64;
        x_re -= (q01_fft[u + (n_uz/2)] * w1_fft[u]) as i64;

        let mut x_im = (q01_fft[u] * w1[u + n_uz/2]) as i64;
        x_im += (q01_fft[u + (n_uz/2)] * w1_fft[u]) as i64;

        let (x_re, z_re) = (x_re.abs(), sign(x_re));
        let (x_im, z_im) = (x_im.abs(), sign(x_im));

        let v = (alpha + q00_fft[u]) as i64;

        if v <= 0 || v >= base_i64.pow(32) || (x_re as i64) >= v*base_i64.pow(32) || (x_im as i64) >= v*base_i64.pow(32)  {
            // fix this also
            return vec![0];
        }

        let y_re = x_re.div_floor(&v);
        let y_im = x_im.div_floor(&v);

        q01_fft[u] = (y_re - 2*z_re*y_re) as i32;
        q01_fft[u + (n_uz/2)] = (y_im - 2*z_im*y_im) as i32;

    }

    let t = ifft(&q01_fft);

    let mut w0: Vec<i32> = vec![0; n_uz];

    for u in 0..n_uz {
        let v = cs0 * h0[u] + t[u];
        let z = (v + cs0) / (2*cs0);

        if z < -((2 as i32).pow(highs0)) || z >= (2 as i32).pow(highs0) {
            return vec![0];
        }

        w0[u] = h0[u] - (2*z);
    }

    // fix this
    return w0;
}

// using this import to floor the quotient after division in the following algorithm
use num::Integer;

pub fn fft(f: &Vec<i32>) -> Vec<i32> {

    let n = f.len();
    let mut f_fft: Vec<i32> = f.clone();
    let mut t = n/2;
    let mut m = 2;

    let mut v0 = 0;

    let f2p31 = (2 as i64).pow(31);
    let f2p32 = (2 as i64).pow(32);

    while m < n {
        v0 = 0;
        for u in 0..m/2{
            let e = delta(u+m);
            let e_re: i64 = e.0 as i64;
            let e_im: i64 = e.1 as i64;

            for v in v0..v0+(t/2) {
                let x1_re: i64 = f_fft[v] as i64;
                let x1_im: i64 = f_fft[v + (n/2)] as i64;

                let x2_re: i64 = f_fft[v + (t/2)] as i64;
                let x2_im: i64 = f_fft[v + (t/2) + (n/2)] as i64;

                let t_re = x2_re*e_re - x2_im*e_im;
                let t_im = x2_re*e_im + x2_im*e_re;

                // f_fft[v] = ((f2p31*x1_re + t_re).div_floor(&f2p32)) as i32;
                f_fft[v] = (((x1_re << 31) + t_re) >> 32) as i32;
                // f_fft[v + (n/2)] = ((f2p31*x1_im + t_im).div_floor(&f2p32)) as i32;
                f_fft[v + (n/2)] = (((x1_im << 31) + t_im) >> 32) as i32;
                // f_fft[v + (t/2)] = ((f2p31*x1_re - t_re).div_floor(&f2p32)) as i32;
                f_fft[v + (t/2)] = (((x1_re << 31) - t_re) >> 32) as i32;
                // f_fft[v + (t/2) + (n/2)] = ((f2p31*x1_im - t_im).div_floor(&f2p32)) as i32;
                f_fft[v + (t/2) + (n/2)] = (((x1_im << 31) - t_im) >> 32) as i32;
            }
            v0 += t;
        }
        t /= 2;
        m *= 2;
    }

    return f_fft;
}

pub fn ifft(f_fft: &Vec<i32>) -> Vec<i32> {

    let n = f_fft.len(); 

    let mut f = f_fft.clone();
    let mut t = 2;
    let mut m = n/2;

    let mut v0 = 0;

    let f2p32 = (2 as i64).pow(32);

    while m > 1 {
        v0 = 0;
        for u in 0..(m/2) {
            let e = delta(u+m);
            let e_re = e.0 as i64;
            let e_im = -(e.0 as i64);

            for v in v0..v0+(t/2) {
                let x1_re = f[v] as i64;
                let x1_im = f[v + (n/2)] as i64;

                let x2_re = f[v + (t/2)] as i64;
                let x2_im = f[v + (t/2) + (n/2)] as i64;

                let t1_re = x1_re + x2_re;
                let t1_im = x1_im + x2_im;

                let t2_re = x1_re - x2_re;
                let t2_im = x1_im - x2_im;

                // f[v] = t1_re.div_floor(&2) as i32;
                f[v] = (t1_re >> 1) as i32;
                // f[v + (n/2)] = t1_im.div_floor(&2) as i32;
                f[v + (n/2)] = (t1_im >> 1) as i32;
                // f[v + (t/2)] = (t2_re * e_re - t2_im*e_im).div_floor(&f2p32) as i32;
                f[v + (t/2)] = ((t2_re * e_re - t2_im*e_im) >> 32) as i32;
                // f[v + (t/2) + (n/2)] = (t2_re*e_im + t2_im*e_re).div_floor(&f2p32) as i32;
                f[v + (t/2) + (n/2)] = ((t2_re*e_im + t2_im*e_re) >> 32 ) as i32;
            }
            v0 += t;
        }
        t *= 2;
        m /=2;
    }

    return f;

}

pub fn polyQnorm(q00: &Vec<i64>, q01: &Vec<i64>, w0: &Vec<i64>, w1: &Vec<i64>, p: i64) -> i64 {
    let n = q00.len();
    let p_u32 = p as u32;

    let q00ntt = ntt(q00.to_vec(), p_u32);
    let q01ntt = ntt(q01.to_vec(), p_u32);
    let w0ntt = ntt(w0.to_vec(), p_u32);
    let w1ntt = ntt(w1.to_vec(), p_u32);

    let mut d: Vec<i64> = Vec::with_capacity(n);
    let mut val = 0;
    for i in 0..n {
        val = modulo(w1ntt[i] * mod_pow(q00ntt[i], p - 2, p), p);
        d.push(val);
    }

    let w1_adj_ntt = ntt(adjoint(&w1), p_u32);
    let mut c: Vec<i64> = Vec::with_capacity(n);
    for i in 0..n {
        val = modulo(d[i] * w1_adj_ntt[i], p);
        c.push(val);
    }

    let mut acc: i64 = c.iter().sum();

    let mut e: Vec<i64> = Vec::with_capacity(n);
    for i in 0..n {
        val = modulo(w0ntt[i] + (d[i] * q01ntt[i]), p);
        e.push(val);
    }

    let e_adj_ntt = nttadj(&e, p_u32);
    for i in 0..n {
        // preprocessing
        let temp1 = modulo(q00ntt[i], p) as u128;
        let temp2 = modulo(e[i], p) as u128;
        let temp3 = modulo(e_adj_ntt[i], p) as u128;
        val = modulo(temp1 * temp2 * temp3, p as u128) as i64;
        acc += val;
    }

    return modulo(acc, p);
}

pub fn poly_div_const_i64(f: &Vec<i64>, c: i64) -> Vec<i64> {
    return f.iter().map(|&x| x / c).collect();
}
pub fn poly_times_const(f: &Vec<i64>, c: i64) -> Vec<i64> {
    return f.iter().map(|&x| x * c).collect();
}

pub fn poly_div_const_f64(f: &Vec<i64>, c: i64) -> Vec<f64> {
    return f.iter().map(|&x| x as f64 / c as f64).collect();
}

fn poly_sub_f64(f: &Vec<f64>, g: &Vec<f64>) -> Vec<f64> {
    let mut h: Vec<f64> = vec![0.0; f.len()];
    for i in 0..h.len() {
        h[i] = f[i] - g[i];
    }
    return h;
}

fn poly_add_f64(f: &Vec<f64>, g: &Vec<f64>) -> Vec<f64> {
    let mut h: Vec<f64> = vec![0.0; f.len()];
    for i in 0..h.len() {
        h[i] = f[i] + g[i];
    }
    return h;
}
