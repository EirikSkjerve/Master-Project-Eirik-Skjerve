use na::*;
use nalgebra as na;

use rand::rngs::StdRng;
use rand::SeedableRng;
use rand::Rng;

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
    m
}

fn nabla_mom4(
    u: &Matrix<f64, Dyn, Dyn, VecStorage<f64, Dyn, Dyn>>,
    w: &Matrix<f64, Dyn, Const<1>, VecStorage<f64, Dyn, Const<1>>>,
) -> Matrix<f64, Const<1>, Dyn, VecStorage<f64, Const<1>, Dyn>> {
    // dot product
    let uw = u * w;
    // power of 3 to each entry
    let uw3 = uw.map(|x| x.powi(3));
    let uw3u = 4.0 * (uw3.clone().transpose() * u) / u.nrows() as f64;
    uw3u
}

pub fn gradient_descent(
    u: Matrix<f64, Dyn, Dyn, VecStorage<f64, Dyn, Dyn>>,
    linv: Matrix<f64, Dyn, Dyn, VecStorage<f64, Dyn, Dyn>>,
    rate: f64,
) -> Matrix<i32, Dyn, Dyn, VecStorage<i32, Dyn, Dyn>> {
    let n = u.ncols();

    // create an empty vec that will keep unique solutions
    let mut solutions: Vec<Matrix<i32, Const<1>, Dyn, VecStorage<i32, Const<1>, Dyn>>> =
        Vec::with_capacity(n);

    let mut iterations = 0;

    let mut rng_seed = rand::thread_rng();
    let mut seed: usize = rng_seed.gen();
    while solutions.len() < n {
        let mut w = gen_u_vec(n, seed);

        seed += 1;
        loop {
            iterations += 1;
            let g = nabla_mom4(&u, &w).transpose();

            let mut w_new = w.clone() - (rate * g);
            let norm_sq: f64 = w_new.iter().map(|&x| x.powi(2)).sum();
            w_new /= norm_sq.sqrt();

            if mom4(&u, &w) <= mom4(&u, &w_new) {
                // round entries
                let vtemp = w.clone().transpose() * linv.clone();
                let v = vtemp.map(|x| x.round() as i32);
                let neg_v = -v.clone();

                if !solutions.contains(&v) && !solutions.contains(&neg_v) {
                    solutions.push(v.clone());
                }
                break;
            } else {
                w = w_new;
            }
        }
    }

    println!(
        "RUST HPP GRADIENT DESCENT COMPLETED IN {} ITERATIONS",
        iterations
    );

    // convert vec of row vectors to matrix object
    let v_mat = Matrix::from_rows(&solutions);

    // return the guess
    v_mat
}
