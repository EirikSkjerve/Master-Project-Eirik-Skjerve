use num::Zero;
use num_complex::{Complex, ComplexFloat};
use std::f64::consts::PI;

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


pub fn sign(x: i32) -> u8 {
    if x < 0{
        return 1;
    }
    return 0;
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

                f_fft[v] = ((f2p31*x1_re + t_re).div_floor(&f2p32)) as i32;
                f_fft[v + (n/2)] = ((f2p31*x1_im + t_im).div_floor(&f2p32)) as i32;
                f_fft[v + (t/2)] = ((f2p31*x1_re - t_re).div_floor(&f2p32)) as i32;
                f_fft[v + (t/2) + (n/2)] = ((f2p31*x1_im - t_im).div_floor(&f2p32)) as i32;

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

                f[v] = t1_re.div_floor(&2) as i32;
                f[v + (n/2)] = t1_im.div_floor(&2) as i32;
                f[v + (t/2)] = (t2_re * e_re - t2_im*e_im).div_floor(&f2p32) as i32;
                f[v + (t/2) + (n/2)] = (t2_re*e_im + t2_im*e_re).div_floor(&f2p32) as i32;
            }
            v0 += t;
        }
        t *= 2;
        m /=2;
    }

    return f;

}
