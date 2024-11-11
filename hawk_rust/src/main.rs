mod fft;
mod fft_constants;
mod ntt;
mod ntt_constants;

mod ntru_solve;
mod rngcontext;
mod utils;

mod parameters;

mod hawk_tests;

mod write_to_file;

mod hawkkeygen;
mod hawksign;
mod hawkverify;
mod verifyutils;

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
