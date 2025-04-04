use hawklib::hawkkeygen::{gen_f_g, hawkkeygen};
use nalgebra::*;

use crate::hawk_sim::hawk_sim_keygen;
use crate::hpp_attack::{measure_res, to_mat};

pub fn compare_keys(n: usize, binv: &DMatrix<i32>) {
    let mut avg_min = 0.0;
    let mut avg_max = 0.0;
    let mut tot_min = f64::INFINITY;
    let mut tot_max = f64::NEG_INFINITY;
    let iterations = 100;
    for _ in 0..iterations {
        let (privkey, _) = hawkkeygen(n, None);
        let (b, binv2) = to_mat(&privkey);

        for i in 0..2 {
            let b1 = binv2.column(i * n + 1).into_owned().map(|x| x as i32);

            let (min, max) = measure_res(&b1, binv);
            avg_min += min as f64 / (2.0 * iterations as f64);
            avg_max += max as f64 / (2.0 * iterations as f64);

            if min < tot_min {
                tot_min = min
            }
            if max > tot_max {
                tot_max = max
            }
        }
    }

    println!("Avg min: {avg_min} \nAvg max: {avg_max}");
    println!("Total min: {tot_min} \nTotal max: {tot_max}");
}

pub fn compare_keys_simulated(n: usize) {
    let ((b, binv), q) = hawk_sim_keygen(n);

    let mut avg_min = 0.0;
    let mut avg_max = 0.0;
    let mut tot_min = f64::INFINITY;
    let mut tot_max = f64::NEG_INFINITY;
    let iterations = 1000;
    for _ in 0..iterations {
        let ((b1, binv1), q1) = hawk_sim_keygen(n);

        for i in 0..2 {
            let b1_col = binv1.column(i * n + 1).into_owned().map(|x| x as i32);

            let (min, max) = measure_res(&b1_col, &binv.map(|x| x as i32));
            avg_min += min as f64 / (2.0 * iterations as f64);
            avg_max += max as f64 / (2.0 * iterations as f64);

            if min < tot_min {
                tot_min = min
            }
            if max > tot_max {
                tot_max = max
            }
        }
    }

    println!("Avg min: {avg_min} \nAvg max: {avg_max}");
    println!("Total min: {tot_min} \nTotal max: {tot_max}");
}
