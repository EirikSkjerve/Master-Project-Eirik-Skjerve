use nalgebra::*;
use rayon::prelude::*;

use rand::distributions::{Distribution, Uniform};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use std::io::{stdout, Write};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use peak_alloc::PeakAlloc;

static PEAK_ALLOC: PeakAlloc = PeakAlloc;

static EPSILON: f32 = 1e-10;
static TOLERANCE: f32 = 1e-4;
static BETA1: f32 = 0.9;
static BETA2: f32 = 0.999;
static DELTA: f32 = 0.01;

fn gradient_optimize(u: &DMatrix<f32>, c: &DMatrix<f32>, descent: bool) -> Option<DVector<f32>> {
    // perform gradient descent if parameter <descent> is set to true
    // otherwise do gradient ascent.

    let n = u.nrows();
    let seed = rand::random::<u64>();
    if descent {
        println!("Running gradient descent");
    } else {
        println!("Running gradient ascent")
    }
    println!("Seed: {seed}");
    let mut rng = StdRng::seed_from_u64(seed);

    // random starting point on the unit circle
    let mut w = get_rand_w(n, &mut rng);
    let mut mom4_w = mom4(&w, &u);

    // initialize variables

    let mut m_t = DVector::<f32>::zeros(n);
    let mut v_t = DVector::<f32>::zeros(n);
    let mut m_est = DVector::<f32>::zeros(n);
    let mut v_est = DVector::<f32>::zeros(n);
    let mut t = 0;

    loop {
        t += 1;
        println!("\nAt iteration {t}");

        // compute gradient of fourth moment
        let g = grad_mom4_par(&w, &u);

        m_t = (BETA1 * &m_t) + (1.0 - BETA1) * &g;
        v_t = (BETA2 * v_t) + (1.0 - BETA2) * &g.map(|x| x.powi(2));
        m_est = &m_t / (1.0 - BETA1.powi(t));
        v_est = &v_t / (1.0 - BETA2.powi(t));

        // based on descent/ascent we need to move in different directions
        let mut wnew = match descent {
            true => &w - (m_est.zip_map(&v_est, |mi, vi| (DELTA * mi) / (vi.sqrt() + EPSILON))),
            false => &w + (m_est.zip_map(&v_est, |mi, vi| (DELTA * mi) / (vi.sqrt() + EPSILON))),
        };

        // normalize new w to put it back on the unit circle
        wnew = wnew.normalize();

        // compute 4th moment of new w
        let mom4_wnew = mom4(&wnew, &u);

        // measure difference between old and new
        let mom4_diff = mom4_wnew - mom4_w;
        let w_diff = (&w - &wnew).norm();

        println!("Mom4(w)    : {mom4_w}");
        println!("Mom4(w_new): {mom4_wnew}");
        println!("Diff       : {mom4_diff}");

        // check requirements for extremum point
        // if mom4_diff.abs() <= TOLERANCE || w_diff <= TOLERANCE || g.norm() <= TOLERANCE {
        //     println!("Possibly at an extremum: rate of change negligible. \n");
        //     return Some(w);
        // }

        // in the case of gradient descent
        if descent && mom4_wnew >= mom4_w {
            println!("Possibly at an extremum: mom4 started increasing. \n");
            return Some(w);
        }

        // in the case of gradient ascent
        if !descent && mom4_wnew <= mom4_w {
            println!("Possibly at an extremum: mom4 started increasing. \n");
            return Some(w);
        }

        // before next iteration, update current w and mom4(w)
        w = wnew;
        mom4_w = mom4_wnew;

        if t >= 100 {
            println!("Too many steps - probably?");
            return Some(w);
        }
        println!("Max memory used so far: {} GB", PEAK_ALLOC.peak_usage_as_gb());
        println!("");
    }
}

pub fn gradient_descent(samples: &DMatrix<f32>, c: &DMatrix<f32>) -> Option<DVector<f32>> {
    if let result = gradient_optimize(samples, c, true) {
        return result;
    } else {
        return None;
    }
}

pub fn gradient_ascent(samples: &DMatrix<f32>, c: &DMatrix<f32>) -> Option<DVector<f32>> {
    if let result = gradient_optimize(samples, c, false) {
        return result;
    } else {
        return None;
    }
}

fn mom4(w: &DVector<f32>, samples: &DMatrix<f32>) -> f32 {
    let dot: DVector<f32> = (samples.transpose() * w).map(|x| x.powi(4));
    dot.mean()
}

fn grad_mom4(w: &DVector<f32>, samples: &DMatrix<f32>) -> DVector<f32> {
    let n = w.nrows();
    let t = samples.ncols();
    let uw3: DVector<f32> = (samples.transpose() * w).map(|x| x.powi(3));
    let uw3u = DMatrix::from_fn(n, t, |i, j| uw3[j] * samples[(i, j)]);
    let g = 4.0 * uw3u.column_mean();
    g
}

fn grad_mom4_par(w: &DVector<f32>, samples: &DMatrix<f32>) -> DVector<f32> {
    let n = w.nrows();
    let t = samples.ncols();

    // Precompute the dot product of each column of `samples` with `w`
    let uw: Vec<f32> = (0..t)
        .into_par_iter() // Parallelize over columns
        .map(|j| samples.column(j).dot(w))
        .collect();

    // Compute uw^3 for each column
    let uw3: Vec<f32> = uw.par_iter().map(|x| x.powi(3)).collect();

    // Compute the gradient in parallel
    let mut g = Arc::new(Mutex::new(DVector::zeros(n)));
    (0..n).into_par_iter().for_each(|i| {
        let mut sum = 0.0;
        for j in 0..t {
            sum += uw3[j] * samples[(i, j)];
        }
        g.lock().unwrap()[i] = 4.0 * sum / (t as f32);
    });

    Arc::try_unwrap(g)
        .expect("Could not unpack gradient g")
        .into_inner()
        .unwrap()
}

fn get_rand_w(n: usize, rng: &mut StdRng) -> DVector<f32> {
    //  outputs some randomly generated w on the unit circle

    // define uniform distribution
    let dist = Uniform::from(-10.0..10.0);

    // initialize empty vector to store the samples
    let mut rnd_bytes: Vec<f32> = Vec::with_capacity(n);

    // sample n times
    for _ in 0..n {
        rnd_bytes.push(dist.sample(rng));
    }

    // load random number into DVector
    let mut w = DVector::from_vec(rnd_bytes);

    // normalize the w vector
    w = w.normalize();
    w
}
