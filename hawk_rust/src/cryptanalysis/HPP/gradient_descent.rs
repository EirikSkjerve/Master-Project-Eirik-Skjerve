use crate::rngcontext::{get_random_bytes, RngContext};
use na::*;
use nalgebra as na;
use std::time::Instant;

use rand::distributions::{Distribution, Uniform};
use rand::rngs::StdRng;
use rand::SeedableRng;

use std::collections::HashSet;

use crate::cryptanalysis::HPP::hpp::get_uni_slice_float;

/// generate some uniformly distributed unit vector with entries in range [-1..1]
fn gen_u_vec(n: usize, seed: usize) -> Matrix<f64, Dyn, Const<1>, VecStorage<f64, Dyn, Const<1>>> {
    // generate some random numbers
    let z = get_uni_slice_float(n, 1, &mut StdRng::seed_from_u64(seed as u64));
    // compute the squared norm and norm
    let z_norm_sq: f64 = z.iter().map(|&x| x.powf(2.0)).sum();
    let z_norm: f64 = (z_norm_sq as f64).sqrt();
    // divide each entry by the norm to obtain a unit vector
    let u_z: Vec<f64> = z.iter().map(|&x| x as f64 / z_norm).collect();
    // return the vector
    DVector::from_row_slice(&u_z)
}

fn mom4(
    u: &Matrix<f64, Dyn, Dyn, VecStorage<f64, Dyn, Dyn>>,
    w: &Matrix<f64, Dyn, Const<1>, VecStorage<f64, Dyn, Const<1>>>,
) -> f64 {

    // dot product
    let uw = u * w;
    // power of 4 to each entry
    let uw4 = uw.map(|x| x.powi(4));
    // the mean value
    let m = uw4.map(|x| x).sum() / (u.nrows() as f64);
    // eprintln!("from mom4 returning m={}", m);
    m
}

fn nabla_mom4(
    u: &Matrix<f64, Dyn, Dyn, VecStorage<f64, Dyn, Dyn>>,
    w: &Matrix<f64, Dyn, Const<1>, VecStorage<f64, Dyn, Const<1>>>,
) -> Matrix<f64, Const<1>, Dyn, VecStorage<f64, Const<1>, Dyn>>{

    // eprintln!("input to nabla_mom4: U: {u} \nw: {w}");
    // println!("Doing nablamom4");
    // dot product
    let uw = u * w;
    // println!("uw shape({}, {})", uw.nrows(), uw.ncols());
    // power of 3 to each entry
    let uw3 = uw.map(|x| x.powi(3));
    let uw3u = 4.0*(uw3.clone().transpose()*u) / u.nrows() as f64;
    // eprintln!("from nabla_mom4 returning uw3u: {uw3u}");
    uw3u

}

pub fn gradient_descent(
    u: Matrix<f64, Dyn, Dyn, VecStorage<f64, Dyn, Dyn>>,
    linv: Matrix<f64, Dyn, Dyn, VecStorage<f64, Dyn, Dyn>>,
    rate: f64,
) {

    eprintln!("U: {u:.2}");
    eprintln!("l⁻¹: {linv:.2}");

    let n = u.ncols();

    // create an empty hash-set that will keep unique solutions
    let mut solutions: HashSet<Matrix<i32, Const<1>, Dyn, VecStorage<i32, Const<1>, Dyn>>> = HashSet::with_capacity(n);

    let mut iterations = 0;

    let mut seed_ctr = 1337;
    while solutions.len() < n {
        let mut w = gen_u_vec(n, seed_ctr);

        // eprintln!("Sample w: {}", w);

        seed_ctr += 1;
        loop{
            iterations += 1;
            let g = nabla_mom4(&u, &w).transpose();

            let mut w_new = w.clone() - (rate*g);
            let temp: f64 = w_new.iter().map(|&x| x.powf(2.0)).sum();
            w_new /= temp.sqrt();

            if mom4(&u, &w) <= mom4(&u, &w_new) {
                // println!("Inside loop");
                // round entries
                // eprintln!("vtemp: {}, {}", linv.nrows(), linv.ncols());
                let vtemp = w.clone().transpose() * linv.clone();
                let v = vtemp.map(|x| x.round() as i32);
                // eprintln!("v: {v}");
                let neg_v = -v.clone();

                if !solutions.contains(&v) && !solutions.contains(&neg_v) {
                    solutions.insert(v.clone());
                }
                break;
            }
            else {
                w = w_new;
            }

        }
    }

    println!("RUST HPP GRADIENT DESCENT COMPLETED IN {} ITERATIONS", iterations);

    let mut retvec = vec![];
    for (i, sol) in solutions.iter().enumerate() {
        retvec.push(sol);
        eprintln!("{:?}", sol);
    }

}
