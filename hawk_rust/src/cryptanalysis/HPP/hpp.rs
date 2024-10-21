use nalgebra as na;
use na::*;
use crate::rngcontext::RngContext;

const NUM_SAMPLES: usize = 10000;
const N: usize = 10;

fn get_uni_slice_int(n: usize, entry_bound: usize) -> Vec<i16> {

    // get some random bytes (i.e. integers in [0..255], that are uniformly distributed
    // using fixed seed
    let mut rng = RngContext::new(&[4]);
    let rnd_bytes_u8 = rng.random(n);

    // make copy with reduced entries mod n, and convert them to i16 for possible negation
    let mut rnd_bytes_i16: Vec<i16> = rnd_bytes_u8.iter().map(|&x| (x % (entry_bound+1) as u8) as i16).collect();
    for (i, rb) in rnd_bytes_u8.iter().enumerate(){
        // negating 50% of the entries
        if *rb >= 128 {
            rnd_bytes_i16[i] = -rnd_bytes_i16[i];
        }
    }
    rnd_bytes_i16

}

fn get_uni_slice_float(n: usize, dist_bound: usize) -> Vec<f64> {
    
    // get some random bytes, uniformly distributed
    // using fixed seed
    let mut rng = RngContext::new(&[1]);
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
fn gen_sec_mat(degree: usize, entry_bound: usize) -> Matrix<i16, Dyn, Dyn, VecStorage<i16, Dyn, Dyn>>{

    // generate slice of uniformly random integers between -entry_bound..entry_bound
    let uni_slice = get_uni_slice_int(degree*degree, entry_bound);

    // create a matrix from this slice
    let sec_v = DMatrix::from_column_slice(degree, degree, &uni_slice);
    
    // return secret matrix v
    sec_v
}


pub fn run_hpp_attack() {

    let entry_bound = 1;
    let dist_bound = 1;

    // create random secret vector v
    let sec_v = gen_sec_mat(N, entry_bound);
    eprintln!("{sec_v}");
    // create a copy of v with floating point numbers
    let sec_v_f = sec_v.map(|x| x as f64); 
    
    // just hardcoding this for now
    let ex2 = 0.333;

    let mut uni_samples = vec![];

    for _ in 0..NUM_SAMPLES {
        let x = get_uni_slice_float(N, dist_bound);
        let x_vec = DVector::from_row_slice(&x);
        let y_vec = x_vec.transpose()*sec_v_f.clone();
        // eprintln!("x: {x_vec:.4}");
        uni_samples.push(y_vec);
    }

    let pub_y = DMatrix::from_rows(&uni_samples);
    
    let g_approx_f = 
        (1.0/ex2)*
        (pub_y.transpose() * pub_y)* 
        (1.0/NUM_SAMPLES as f64);

    // eprintln!("{g_approx_f:.2}");

    let g_approx = g_approx_f.map(|x| x.round());
    // let g_approx = g_approx_f.clone();

    let mut g_approx_inverse = g_approx.clone();
    match g_approx_inverse.clone().try_inverse(){
        Some(g_inv) => {
            g_approx_inverse = g_inv.clone();
            // use this for good formatting
            eprintln!("g_approx: {}",g_approx);
            eprintln!("g_approx_inv: {g_approx_inverse:.4}");
            let id_mat = g_approx * g_approx_inverse.clone();
            eprintln!("{id_mat:.2}");

        }
        None => {
            println!("Not invertible");
        }
    }

    let l: Cholesky<f64, Dyn>;
    if let Some(cholesky) = g_approx_inverse.clone().cholesky() {
        l = cholesky
    } else {
        println!("Could not decompose matrix");
        return
    };
    
    let linv = l.clone().inverse();

}
