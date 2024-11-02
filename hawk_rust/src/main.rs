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

mod linalg;
mod poly;

mod cryptanalysis;
mod write_to_file;

use cryptanalysis::HPP::hpp::run_hpp_attack;
use hawk_tests::test_all;

// memory measurement
use peak_alloc::PeakAlloc;

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;


fn main() {

    // run several instances of hawk 256, 512, and 1024
    // measures avg. time usage
    test_all();
    // run_hpp_attack();
}
