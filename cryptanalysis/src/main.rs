#![allow(warnings)] // Silences all warnings in the crate
#[macro_use]
extern crate prettytable;
mod collect_signatures;
mod dgd_estimation;
mod dgd_estimation_normalized;
mod hpp;
mod file_utils;

use collect_signatures::{covariance_matrix_estimation, collect_signatures};
use dgd_estimation::{estimate_mem, estimate_mem_all};
use dgd_estimation_normalized::estimate_mem_norm_all;
use hpp::hpp::run_hpp_attack;

use prettytable::{Cell, Row, Table};
use peak_alloc::PeakAlloc;

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

fn main() {
    // run_hpp_attack();
    collect_signatures(7500000, 256);
    covariance_matrix_estimation(7500000,256);
    // estimate_mem_all(500000, true);
    // estimate_mem_norm_all(3000000, false);
    println!("Max memory usage total in this run: {} GB", PEAK_ALLOC.peak_usage_as_gb());
}
