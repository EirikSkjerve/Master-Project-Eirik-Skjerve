use crate::rngcontext::{get_random_bytes, shake256x4, RngContext};
use hawklib::ntru_solve::ntrusolve;
use hawklib::utils::rot_key;
use hawklib::hawksign::sample;

use rand::{thread_rng, Rng};
use rand::seq::SliceRandom;
use rand::distributions::Uniform;
use nalgebra::*;
use rand_distr::{Bernoulli, Distribution};

pub fn hawk_sim_keygen(n: usize) -> ((DMatrix<i64>, DMatrix<i64>), DMatrix<i64>){

    loop {

        let rndbytes = get_random_bytes(20);
        let (f, g) = gen_f_g(n);
        if let Some((bigf, bigg)) = ntrusolve(&f, &g){
            let (b, binv) = to_mat_priv(&f, &g, &bigf, &bigg);
            let q = &b.transpose() * &b;
            return ((b, binv), q);
        }
    }
}

fn random_uniform_vector(n: usize) -> DVector<i64> {
    let mut rng = rand::thread_rng();
    let values: Vec<i64> = (0..n).map(|_| rng.gen_range(-10..10)).collect();
    DVector::from_vec(values)
}
pub fn hawk_sim_sign(n: usize, binv: &DMatrix<i64>) -> DVector<i64> {
    let seed = get_random_bytes(40);
    let mut t = get_random_bytes(2*n);
    t = t.iter().map(|&x| x % 2).collect();
    let x = DVector::<i64>::from_vec(sample(&seed, t, n));

    // as a test, use a uniform distribution here instead 
    let x = random_uniform_vector(2*n);
    let w = binv * x;
    w
}

pub fn hawk_sim_verify(n: usize, sig: &DVector<i64>, q: &DMatrix<i64>) {
    let qnorm = (*(((1/2*n as i64) * sig.transpose() * q * sig).get((0,0))).unwrap() as f64).sqrt();
}

fn to_mat_priv(f: &Vec<i64>, g: &Vec<i64>, bigf: &Vec<i64>, bigg: &Vec<i64>) -> (DMatrix<i64>, DMatrix<i64>){

    let n = f.len();

    let b = rot_key(&f, &g, &bigf, &bigg);
    let binv = rot_key(
        &bigg,
        &g.clone().iter().map(|&x| -x).collect(),
        &bigf.clone().iter().map(|&x| -x).collect(),
        &f,
    );

    // convert them to Nalgebra matrices
    let flatb: Vec<i64> = b.into_iter().flatten().collect();
    let flatbinv: Vec<i64> = binv.into_iter().flatten().collect();

    let b = DMatrix::from_row_slice(2 * n, 2 * n, &flatb);
    let binv = DMatrix::from_row_slice(2 * n, 2 * n, &flatbinv);
    (b, binv)
}

fn gen_f_g(n: usize) -> (Vec<i64>, Vec<i64>){
    let (mut f, mut g) = (Vec::<i64>::new(), Vec::<i64>::new());

    let eta = 2;

    let unif = Uniform::new(0, 1);
    let dist = Bernoulli::new(0.5).unwrap();
    let mut rng = thread_rng();

    for _ in 0..n {

        let mut sample = 0;
        for i in 0..2*eta {
            if dist.sample(&mut rng) {
                sample += 1;
            }
        }
        f.push(sample - eta);
    }
    for _ in 0..n {

        let mut sample = 0;
        for i in 0..2*eta {
            if dist.sample(&mut rng) {
                sample += 1;
            }
        }
        g.push(sample - eta);
    }

    (f, g)
}
