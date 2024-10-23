use nalgebra as na;
use na::*;
use crate::rngcontext::{RngContext, get_random_bytes};

use rand::distributions::{Distribution, Uniform};
use rand::SeedableRng;
use rand::rngs::StdRng;


const NUM_SAMPLES: usize = 5000;
const N: usize = 16;

fn get_uni_slice_int(n: usize, entry_bound: usize, seed: usize) -> Vec<i8> {

    let bound = Uniform::from(-(entry_bound as i8)..(entry_bound as i8));
    let mut rng = StdRng::seed_from_u64(seed as u64);
    let mut rnd_bytes: Vec<i8> = Vec::with_capacity(n);

    for _ in 0..n {
        rnd_bytes.push(bound.sample(&mut rng));
    }

    rnd_bytes

}

fn get_uni_slice_float(n: usize, dist_bound: usize, rng: &mut RngContext) -> Vec<f64> {
    
    // get some random bytes, uniformly distributed
    // using fixed seed
    let rnd_bytes_u8 = rng.random(n);

    // make a copy with normalised entries
    let mut rnd_bytes_f64: Vec<f64> = rnd_bytes_u8
        .iter()
        .map(|&x| (x as f64) / (255.0 / dist_bound as f64))
        .collect();

    for (i, rb) in rnd_bytes_u8.iter().enumerate(){
        // negating 50% of the entries
        if *rb >= 128 {
            rnd_bytes_f64[i] = -rnd_bytes_f64[i];
        }
    }
    rnd_bytes_f64
}

/// generate a secret matrix V than we will try and retrieve
fn gen_sec_mat(degree: usize, entry_bound: usize) -> Matrix<f64, Dyn, Dyn, VecStorage<f64, Dyn, Dyn>>{

    loop {
        // generate slice of uniformly random integers between -entry_bound..entry_bound
        let uni_slice = get_uni_slice_int(degree*degree, entry_bound, 42);
        let uni_slice_f: Vec<f64> = uni_slice.iter().map(|&x| x as f64).collect();

        // create a matrix from this slice
        let sec_v = DMatrix::from_column_slice(degree, degree, &uni_slice_f);

        // perform singular value decomposition to calculate rank of matrix
        let svd_v = sec_v.clone().svd(true, true);
        let rank = svd_v.rank(1e-10);
        
        if rank < degree {
            continue;
        }

        // return secret matrix v
        eprintln!("Returning: {sec_v}");
        return sec_v;
    }
}


pub fn run_hpp_attack() {

    let entry_bound = 1;
    let dist_bound = 1;


    // create random secret vector v
    // let sec_v = gen_sec_mat(N, entry_bound);

    // using a fixed matrix V
    // let v_data: [i16; 256] = [1, -1, -1, -1, -1,  0,  0,  1,  1,  1,  0, -1, -1,  0,  1,  0, 0,  0, -1,  0, -1,  1, -1,  1,  1,  0, -1, -1,  0,  0,  0,  1, 1,  0,  0,  0,  1,  0,  1, -1,  0,  0, -1,  0,  1, -1, -1,  1, 0,  1, -1, -1, -1,  1,  1, -1,  0,  1, -1,  0, -1, -1, -1, -1, 0,  1, -1,  1,  1,  1,  0,  0,  1,  1,  0,  1, -1, -1, -1, -1, 0, -1, -1, -1,  0, -1, -1, -1,  0,  1, -1, -1,  0,  0, -1,  0, 0, -1,  0,  1,  0,  1,  1,  0,  1, -1,  1,  0,  0,  1,  1,  1, -1, -1,  1, -1, -1,  1,  0,  1, -1,  0,  1,  1,  0,  0,  1, -1, -1,  1,  0,  1,  1,  1,  0,  0,  0,  1, -1,  1,  1,  1,  0, -1, 0,  1,  1, -1, -1, -1,  1,  1,  0,  1,  1,  0,  1, -1,  1,  1, -1, -1,  1, -1,  0,  0,  1,  1, -1,  1,  1,  0,  1,  0,  1,  1, 0, -1,  0,  1,  1,  1, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0, 1,  1,  1,  0,  1, -1,  0, -1,  0,  0, -1,  1,  0,  0,  0, -1, 0,  1,  1,  1,  1,  1,  1,  1, -1,  0,  1,  0, -1,  1,  0,  0, 0,  1,  1,  0,  1,  0,  1, -1,  0, -1,  1,  0,  1,  1, -1, -1, 1,  1,  1,  0,  1, -1, -1,  1,  0,  1, -1,  0,  0,  1,  1,  1];

    let sec_v_f = gen_sec_mat(N, entry_bound);

    loop {

        let init_seed = get_random_bytes(5);
        let init_seed = [213, 66, 229, 146, 87];
        let mut rng = RngContext::new(&init_seed);

        // just hardcoding this for now
        let ex2 = 0.333;

        // empty vector storing samples
        let mut uni_samples = vec![];


        // generate a bunch of samples (that are uniformly distributed)
        // and multiply them with secret matrix v
        for i in 0..NUM_SAMPLES {
            // this is slow as fuck, but oh well
            let x = get_uni_slice_float(N, dist_bound, &mut rng);
            let x_vec = DVector::from_row_slice(&x);
            let y_vec = x_vec.transpose() * sec_v_f.clone();
            uni_samples.push(y_vec);
        }

        let pub_y = DMatrix::from_rows(&uni_samples);
        
        // approximation of Gram Matrix
        let g_approx_f = 
            (1.0/ex2)*
            (pub_y.transpose() * pub_y.clone())* 
            (1.0/NUM_SAMPLES as f64);


        let g_approx = g_approx_f.map(|x| x.round());
        // eprintln!("g_approx: {g_approx}");
        // let g_approx = g_approx_f.clone();

        let mut g_approx_inverse = g_approx.clone();
        match g_approx_inverse.clone().try_inverse(){
            Some(g_inv) => {
                g_approx_inverse = g_inv;
                // use this for good formatting
                // eprintln!("g_approx: {}",g_approx);
                // eprintln!("g_approx_inv: {g_approx_inverse:.4}");
                // let id_mat = g_approx * g_approx_inverse.clone(); //.map(|x| x.abs() as i8)
                // eprintln!("this should be identity matrix {id_mat:.3}");

            }
            None => {
                println!("Not invertible");
            }
        }

        let l: Cholesky<f64, Dyn>;
        if let Some(cholesky) = g_approx_inverse.clone().cholesky() {
            l = cholesky;
            let linv = l.clone().inverse();

            eprintln!("L^(-1): {linv:.3}");
        } else {
            println!("Could not decompose matrix");
            continue
        };

        let u = pub_y*l.l()/dist_bound as f64;
        // eprintln!("{u:.5}");
        eprintln!("init_seed: {:?}", init_seed);

        return;
        
    }

}
