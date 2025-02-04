use nalgebra::*;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rand::distributions::{Distribution, Uniform};

use std::io::{stdout, Write};

use peak_alloc::PeakAlloc;

static PEAK_ALLOC: PeakAlloc = PeakAlloc;

static EPSILON: f64 = 1e-10;
static BETA1: f64 = 0.9;
static BETA2: f64 = 0.999;
static DELTA: f64 = 0.7;

fn gradient_optimize(u: &DMatrix<f64>, c: &DMatrix<f64>, descent: bool) -> Option<DVector<f64>> {
    // perform gradient descent if parameter <descent> is set to true
    // otherwise do gradient ascent.

    let n = u.nrows();
    let mut rng = StdRng::seed_from_u64(rand::random::<u64>());

    // random starting point on the unit circle
    let mut w = get_rand_w(n, &mut rng);
    let mut mom4_w = mom4(&w, &u);

    // initialize variables

    let mut m_t = DVector::<f64>::zeros(n);
    let mut v_t = DVector::<f64>::zeros(n);
    let mut t = 0;

    loop {
        t += 1;
        let g = grad_mom4(&w, &u); // get gradient

        m_t = (BETA1*&m_t) + (1.0-BETA1)*&g;
        v_t = (BETA2*v_t) + (1.0-BETA2)*&g.map(|x| x.powi(2));
        let m_est = &m_t / (1.0 - BETA1.powi(t));
        let v_est = &v_t / (1.0 - BETA2.powi(t));

        // based on descent/ascent we need to move in different directions
        let mut wnew = match descent {
            true => &w - (m_est.zip_map(&v_est, |mi, vi| (DELTA*mi)/(vi.sqrt() + EPSILON))),
            false => &w + (m_est.zip_map(&v_est, |mi, vi| (DELTA*mi)/(vi.sqrt() + EPSILON))),
        };
        
        // normalize new w to put it back on the unit circle
        wnew = wnew.normalize();

        // compute 4th moment of new w
        let mom4_wnew = mom4(&wnew, &u);

        // measure difference between old and new
        let mom4_diff = (mom4_w - mom4_wnew).abs();
        let w_diff = (&w - &wnew).norm();

        // check requirements for extremum point
        if mom4_diff <= EPSILON || w_diff <= EPSILON || g.norm() <= EPSILON {
            println!("Possibly at an extremum here...");
            return Some(w);
        }

        // before next iteration, update current w and mom4(w)
        w = wnew;
        mom4_w = mom4_wnew;
    }

}

pub fn gradient_descent(samples: &DMatrix<f64>, c: &DMatrix<f64>) -> Option<DVector<f64>> {
    if let result = gradient_optimize(samples, c, true) {
        return result;
    } else {return None;}
}

pub fn gradient_ascent(samples: &DMatrix<f64>, c: &DMatrix<f64>) -> Option<DVector<f64>> {
    if let result = gradient_optimize(samples, c, false) {
        return result;
    } else {return None;}
}

fn get_rand_w(n: usize, rng: &mut StdRng) -> DVector<f64> {
    //  outputs some randomly generated w on the unit circle

    // define uniform distribution
    let dist = Uniform::from(-10.0..10.0);

    // initialize empty vector to store the samples
    let mut rnd_bytes: Vec<f64> = Vec::with_capacity(n);

    // sample n times
    for _ in 0..n {
        rnd_bytes.push(dist.sample(rng));
    }

    // load random number into DVector
    let mut w = DVector::from_vec(rnd_bytes);

    // normalize the w vector
    w /= w.norm();
    w
}

fn mom4(w: &DVector<f64>, samples: &DMatrix<f64>) -> f64 {

    // naive method
    // let n = samples.nrows();
    // let t = samples.ncols();
    // let mut res: f64 = 0.0;
    // for i in (0..t) {
    //     let dot = samples.column(i).dot(w).powi(4);
    //     // println!("dot: {dot}");
    //     // println!("\nnorm s: {}", samples.column(i).norm());
    //     // println!("norm w: {}", w.norm());
    //     res += dot;
    // }
    // res /= t as f64;
    // println!("Res: {res}");
    // res
    // estimate 4th moment given samples and vector w
    // compute <u,w>^4
    let dot: DVector<f64> = (samples.transpose() * w).map(|x| x.powi(4));
    // // println!("Res 2: {}", dot.mean());
    // // compute mean of above, and return result
    dot.mean()
}

fn grad_mom4(w: &DVector<f64>, samples: &DMatrix<f64>) -> DVector<f64> {
    let n = w.nrows();
    let t = samples.ncols();
    // estimate gradient of 4th moment given samples and vector w
    // compute 4(<u, w>^3 * u)
    // println!("Shape of U: ({}, {})", samples.nrows(), samples.ncols());
    // println!("Shape of w: ({}, {})", w.nrows(), w.ncols());
    // dot product
    // power of 3 to each entry
    // println!("Shape of uw: ({}, {})", uw3.nrows(), uw3.ncols());
    let uw3: DVector<f64> = (samples.transpose() * w).map(|x| x.powi(3));
    // let uw3u: DVector<f64> = (4.0 * (samples * uw3) / samples.nrows() as f64);
    // println!("Shape of uw3u: ({}, {})", uw3u.nrows(), uw3u.ncols());
    // let g: DVector<f64> = (4.0 * (samples * uw3) / samples.ncols() as f64);
    let uw3u = DMatrix::from_fn(n, t, |i, j| uw3[j] * samples[(i, j)]);
    let g = 4.0*uw3u.column_mean();
    // println!("Shape of g: ({}, {})", g.nrows(), g.ncols());
    g
}
