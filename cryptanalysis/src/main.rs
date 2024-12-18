#![allow(warnings)] // Silences all warnings in the crate
mod collect_signatures;
mod dgd_estimation;
mod hpp;


use collect_signatures::generate_t_signatures;
use dgd_estimation::{estimate_sigma_mem, estimate_sigma_mem_all};
use hpp::hpp::run_hpp_attack;

fn main() {
    // run_hpp_attack();
    // generate_t_signatures(90000,256);
    estimate_sigma_mem_all(1000);
}
