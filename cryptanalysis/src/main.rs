#![allow(warnings)] // Silences all warnings in the crate
#[macro_use]
extern crate prettytable;
mod collect_signatures;
mod dgd_estimation;
mod dgd_estimation_normalized;
mod file_utils;
mod gradient_search;
mod hpp;
mod hpp_attack;
mod hpp_attack_yu;
mod hpp_attack_online;
mod hpp_simulation;
mod procrustes_attack;
mod test_candidate_vec;
mod measure_uy;

mod hawk_sim;
mod rngcontext;

mod compare_keys;

use collect_signatures::{
    collect_signatures_par, covariance_matrix_estimation,
};
use dgd_estimation::{estimate_mem, estimate_mem_all};
use dgd_estimation_normalized::estimate_mem_norm_all;

use hpp_attack::run_hpp_attack;
use hpp_attack_yu:: run_hpp_attack_yu;
use hpp_simulation::run_hpp_sim;
use hpp_attack_online::hpp_attack_online;

use procrustes_attack::*;

use peak_alloc::PeakAlloc;
use prettytable::{Cell, Row, Table};

use std::time::{Duration, Instant};
use std::env;

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

fn main() {
    let args: Vec<String> = env::args().collect();
    let mode: &str = &args[1];
    let t: usize = args[2]
        .parse()
        .expect("Invalid input for number of samples");
    let n: usize = args[3].parse().expect("Invalid input for Hawk degree");

    let start = Instant::now();
    match mode { 
        "measure" => estimate_mem_norm_all(t, true),
        "collect" => {collect_signatures_par(t, n, true);},
        "attack" => run_hpp_attack(t, n),
        "attack_yu" => run_hpp_attack_yu(t, n),
        "attack_o" => hpp_attack_online(t, n),
        "sim" => run_hpp_sim(t, n),
        _ => println!("No action matches your input"),
    };

    println!(
        "Max memory usage total in this run: {} GB",
        PEAK_ALLOC.peak_usage_as_gb()
    );

    println!("
             Time used: {:?}",
             start.elapsed());
}
