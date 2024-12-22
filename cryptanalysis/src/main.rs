#![allow(warnings)] // Silences all warnings in the crate
#[macro_use]
extern crate prettytable;
mod collect_signatures;
mod dgd_estimation;
mod hpp;

use collect_signatures::generate_t_signatures;
use dgd_estimation::{estimate_mem, estimate_mem_all};
use hpp::hpp::run_hpp_attack;

use prettytable::{Cell, Row, Table};
use peak_alloc::PeakAlloc;

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

fn main() {
    // run_hpp_attack();
    generate_t_signatures(1000000,256);
    // estimate_mem_all(500000, true);
    //let current_mem = PEAK_ALLOC.current_usage_as_mb();
	// println!("This program currently uses {} MB of RAM.", current_mem);
	let peak_mem = PEAK_ALLOC.peak_usage_as_gb();
	println!("The max amount that was used {} gb", peak_mem);
}
