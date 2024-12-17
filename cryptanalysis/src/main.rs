#![allow(warnings)] // Silences all warnings in the crate
mod collect_signatures;
mod dgd_estimation;
mod hpp;


use collect_signatures::generate_t_signatures;
use dgd_estimation::{estimate_sigma_mem, estimate_sigma_time};
use hpp::hpp::run_hpp_attack;

fn main() {
    // run_hpp_attack();
    // generate_t_signatures(900000,256);
    estimate_sigma_mem(10000000, 256);
}
