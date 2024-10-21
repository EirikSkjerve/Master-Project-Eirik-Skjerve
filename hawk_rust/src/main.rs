mod fft;
mod fft_constants;
mod ntt;
mod ntt_constants;

mod hawk_tests;

mod ntru_solve;
mod rngcontext;
mod utils;

mod parameters;

mod hawk1024;
mod hawk256;
mod hawk512;

mod compression;

mod poly;
mod linalg;

mod cryptanalysis;
mod write_to_file;

use hawk_tests::test_all;
use cryptanalysis::HPP::hpp::run_hpp_attack;

// memory measurement
use peak_alloc::PeakAlloc;

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

/*
   Driver code for HAWK implementation.

   Integer types:
   n: 256-1024, use u16,
   log n: 8-10, use u8,
   vectors of bits: 0/1, use u8,
   vectors f and g: initially low numbers, but is changed fast.
       set as Vec<i64>, will have negative values

*/

fn main() {
    // test_all();
    run_hpp_attack();
}
