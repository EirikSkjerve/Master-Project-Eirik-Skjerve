#![allow(warnings)] // Silences all warnings in the crate
#[macro_use]
extern crate prettytable;
mod collect_signatures;
mod dgd_estimation;
mod dgd_estimation_normalized;
mod file_utils;
mod hpp;
mod hpp_attack;

use collect_signatures::{collect_signatures, covariance_matrix_estimation};
use dgd_estimation::{estimate_mem, estimate_mem_all};
use dgd_estimation_normalized::estimate_mem_norm_all;
// use hpp::hpp::run_hpp_attack;
use hpp_attack::run_hpp_attack;

use peak_alloc::PeakAlloc;
use prettytable::{Cell, Row, Table};

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

fn main() {
    let t = 3000000;
    let n = 256;
    // run_hpp_attack();
    // covariance_matrix_estimation(7500000,256);
    // estimate_mem_all(500000, true);
    // estimate_mem_norm_all(40000, true);

    // collect_signatures(t, n);
    run_hpp_attack(t, n);

    println!(
        "Max memory usage total in this run: {} GB",
        PEAK_ALLOC.peak_usage_as_gb()
    );
}
