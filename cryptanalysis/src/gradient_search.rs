use nalgebra::*;
use rayon::prelude::*;

use rand_distr::{Distribution, Normal, Uniform};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use std::io::{stdout, Write};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use peak_alloc::PeakAlloc;

static PEAK_ALLOC: PeakAlloc = PeakAlloc;

static EPSILON: f64 = 1e-10;
// static TOLERANCE: f64 = 1e-4;
static BETA1: f64 = 0.9;
static BETA2: f64 = 0.999;
static DELTA: f64 = 0.7;
static VANILLA_DELTA: f64 = 0.1;

pub fn gradient_descent_vanilla(u: &DMatrix<f64>) -> Option<DVector<f64>> {
    let n = u.nrows();
    let seed = rand::random::<u64>();

    let mut rng = StdRng::seed_from_u64(seed);

    let mut w = get_rand_w(n, &mut rng);

    let mut mom4_w = mom4_par(&w, &u);

    let mut t = 0;
    loop {
        t += 1;
        println!("Iteration {t}");
        let g = grad_mom4_par(&w, &u);
        let wnew = &w - (VANILLA_DELTA * g);
        let wnew = wnew.normalize();

        let mom4_wnew = mom4_par(&wnew, &u);
        println!("Mom4(w):    {mom4_w}");

        if mom4_wnew >= mom4_w {
            println!("Mom4 started increasing. Returning");
            return Some(w);
        }
        w = wnew;
        mom4_w = mom4_wnew;
    }
    return None;
}

pub fn gradient_ascent_vanilla(u: &DMatrix<f64>, delta: f64) -> Option<DVector<f64>> {
    let n = u.nrows();
    let seed = rand::random::<u64>();

    let mut rng = StdRng::seed_from_u64(seed);

    let mut w = get_rand_w(n, &mut rng);

    let mut mom4_w = mom4_par(&w, &u);

    let mut t = 0;
    let mut inner_retries = 0;
    loop {
        t += 1;
        println!("Iteration {t}");
        let g = grad_mom4_par(&w, &u);
        let wnew = &w + (VANILLA_DELTA * g);
        let wnew = wnew.normalize();

        let mom4_wnew = mom4_par(&wnew, &u);
        println!("Mom4(w):    {mom4_w}");
        println!("Mom4(wnew): {mom4_wnew}");

        if mom4_wnew <= mom4_w {
            println!("Mom4 started decreasing. Returning");
            return Some(w);
        }
        w = wnew;
        mom4_w = mom4_wnew;
    }
    return None;
}

fn gradient_optimize(u: &DMatrix<f64>, descent: bool, solution: Option<&DVector<f64>>) -> Option<DVector<f64>> {
    // perform gradient descent if parameter <descent> is set to true
    // otherwise do gradient ascent.

    let n = u.nrows();
    let seed = rand::random::<u64>();

    println!("Seed: {seed}");
    let mut rng = StdRng::seed_from_u64(seed);

    // random starting point on the unit circle
    let mut w = get_rand_w(n, &mut rng);

    // if a solution vector is provided, set w to this
    // for testing purposes
    if solution.is_some() {
        w = solution.unwrap().clone();
    }
    let mut mom4_w = mom4_par(&w, &u);

    // initialize variables

    let mut m_t = DVector::<f64>::zeros(n);
    let mut v_t = DVector::<f64>::zeros(n);
    let mut m_est = DVector::<f64>::zeros(n);
    let mut v_est = DVector::<f64>::zeros(n);
    let mut t = 0;

    loop {
        t += 1;
        // println!("\nAt iteration {t}");

        // compute gradient of fourth moment
        let g = grad_mom4_par(&w, &u);
        let grad_norm = g.norm();
        // println!("Norm of gradient g: {grad_norm}");

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
        let mom4_wnew = mom4_par(&wnew, &u);

        // measure difference between old and new
        let mom4_diff = mom4_wnew - mom4_w;
        let w_diff = (&w - &wnew).norm();

        println!("Mom4(w)    : {mom4_w}");
        println!("Mom4(w_new): {mom4_wnew}");
        // println!("Diff       : {mom4_diff}");

        // in the case of gradient descent
        if descent && mom4_wnew > mom4_w {
            println!("Possibly at an extremum: mom4 started increasing. \n");
            // return Some(w.map(|x| x.round()));
            return Some(w);
        }

        // in the case of gradient ascent
        if !descent && mom4_wnew < mom4_w {
            println!("Possibly at an extremum: mom4 started decreasing. \n");
            // return Some(w.map(|x| x.round()));
            return Some(w);
        }

        // before next iteration, update current w and mom4(w)
        w = wnew;
        mom4_w = mom4_wnew;

        // if t >= 100 {
        //     println!("Too many steps - probably?");
        //     return Some(w);
        // }
        // println!("Max memory used so far: {} GB", PEAK_ALLOC.peak_usage_as_gb());
        // println!("");
    }
}

pub fn gradient_descent(samples: &DMatrix<f64>, solution: Option<&DVector<f64>>) -> Option<DVector<f64>> {
    if let result = gradient_optimize(samples, true, solution) {
        return result;
    }
    None
}

pub fn gradient_ascent(samples: &DMatrix<f64>, solution: Option<&DVector<f64>>) -> Option<DVector<f64>> {
    if let result = gradient_optimize(samples, false, solution) {
        return result;
    }
    None
}

fn mom4(w: &DVector<f64>, samples: &DMatrix<f64>) -> f64 {
    let dot: DVector<f64> = (samples.transpose() * w).map(|x| x.powi(4));
    dot.mean()
}

fn mom4_par(w: &DVector<f64>, samples: &DMatrix<f64>) -> f64 {
    let uw3: Vec<f64> = (0..samples.ncols())
        .into_par_iter()
        .map(|j| samples.column(j).dot(w).powi(4))
        .collect();

    DVector::<f64>::from_vec(uw3).mean()

}

fn grad_mom4(w: &DVector<f64>, samples: &DMatrix<f64>) -> DVector<f64> {
    let n = w.nrows();
    let t = samples.ncols();
    let uw3: DVector<f64> = (samples.transpose() * w).map(|x| x.powi(3));
    let uw3u = DMatrix::from_fn(n, t, |i, j| uw3[j] * samples[(i, j)]);
    let g = 4.0 * uw3u.column_mean();
    g
}

fn grad_mom4_par(w: &DVector<f64>, samples: &DMatrix<f64>) -> DVector<f64> {
    let n = w.nrows();
    let t = samples.ncols();

    // Precompute the dot product of each column of `samples` with `w`
    let uw3: Vec<f64> = (0..t)
        .into_par_iter() // Parallelize over columns
        .map(|j| samples.column(j).dot(w).powi(3))
        .collect();

    // Compute the gradient in parallel
    let mut g = Arc::new(Mutex::new(DVector::zeros(n)));
    (0..n).into_par_iter().for_each(|i| {
        let mut sum = 0.0;
        for j in 0..t {
            sum += uw3[j] * samples[(i, j)];
        }
        g.lock().unwrap()[i] = 4.0 * sum / t as f64;
    });

    // println!("Max memory usage: {} GB", PEAK_ALLOC.peak_usage_as_gb());
    Arc::try_unwrap(g)
        .expect("Could not unpack gradient g")
        .into_inner()
        .unwrap()
}

fn get_rand_w(n: usize, rng: &mut StdRng) -> DVector<f64> {
    //  outputs some randomly generated w on the unit circle

    let sigma = match n / 2 {
        256 => 2.0 * 1.001,
        512 => 2.0 * 1.278,
        _ => 2.0 * 1.299,
    };
    // define uniform distribution
    let dist1 = Normal::new(0.0, 5.0).unwrap();
    let dist2 = Normal::new(0.0, sigma).unwrap();
    let dist3 = Uniform::from(-1.0..1.0);

    // initialize empty vector to store the samples
    let mut rnd_bytes: Vec<f64> = Vec::with_capacity(n);

    // sample n times
    // n/2 for each distribution
    for _ in 0..n/2 {
        rnd_bytes.push(dist3.sample(rng));
    }
    for _ in 0..n/2 {
        rnd_bytes.push(dist3.sample(rng));
    }

    // load random number into DVector
    let mut w = DVector::from_vec(rnd_bytes);
    // eprintln!("w before normalizing: {w}");
    // normalize the w vector
    w = w.normalize();
    // eprintln!("w after normalizing: {w}");
    w
}
