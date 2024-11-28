mod rngcontext;
mod hawk_tests;
use hawk_tests::test_all;

// memory measurement
use peak_alloc::PeakAlloc;

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

fn main() {
    // run several instances of hawk 256, 512, and 1024
    // measures avg. time usage
    test_all();
}
