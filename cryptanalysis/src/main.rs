#![allow(warnings)] // Silences all warnings in the crate
#[macro_use]
extern crate prettytable;
mod collect_signatures;
mod dgd_estimation;
mod dgd_estimation_normalized;
mod file_utils;
mod hpp;
mod hpp_attack;
mod hpp_simulation;
mod test_candidate_vec;

use collect_signatures::{
    collect_signatures, collect_signatures_par, covariance_matrix_estimation,
};
use dgd_estimation::{estimate_mem, estimate_mem_all};
use dgd_estimation_normalized::estimate_mem_norm_all;

use hpp_attack::run_hpp_attack;
use hpp_simulation::run_hpp_sim;

use peak_alloc::PeakAlloc;
use prettytable::{Cell, Row, Table};

use std::env;

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

fn main() {
    let args: Vec<String> = env::args().collect();
    let t: usize = args[1]
        .parse()
        .expect("Invalid input for number of samples");
    let n: usize = args[2]
        .parse()
        .expect("Invalid input for Hawk degree");

    // covariance_matrix_estimation(t, n);
    // estimate_mem_all(500000, true);
    // estimate_mem_norm_all(t, false);

    // collect_signatures_par(t, n);
    run_hpp_attack(t, n);
    // run_hpp_sim(t, n);

    println!(
        "Max memory usage total in this run: {} GB",
        PEAK_ALLOC.peak_usage_as_gb()
    );
}
